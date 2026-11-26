# Beta-Binomial GLM with Cox-Reid Empirical Bayes dispersion shrinkage,
# arbitrary covariate design formulas, likelihood ratio tests, and BH FDR.

const _loggamma = Distributions.SpecialFunctions.loggamma
const _digamma  = Distributions.SpecialFunctions.digamma
const _trigamma = Distributions.SpecialFunctions.trigamma

const _GLM_WHITELIST_EVENTS = Set{String}(["TS", "TE", "AF", "AL", "RI", "CE", "AA", "AD"])

const _THETA_MIN = -4.0   # M = exp(-4)  ≈ 0.0183  (phi ≈ 0.982)
const _THETA_MAX = 11.5   # M = exp(11.5) ≈ 98715  (phi ≈ 1.0e-5, effectively Binomial)
const _MU_EPS    = 1e-9

@inline function _logistic(x::Float64)::Float64
   if x >= 0.0
      return 1.0 / (1.0 + exp(-x))
   else
      ex = exp(x)
      return ex / (1.0 + ex)
   end
end

@inline function _logit(p::Float64)::Float64
   pc = clamp(p, _MU_EPS, 1.0 - _MU_EPS)
   return log(pc / (1.0 - pc))
end

# Continuous Beta-Binomial log-likelihood (omitting constant log-binomial coefficient
# that cancels identically across all beta/M derivatives and LRTs).
function bb_loglik(k::AbstractVector{Float64},
                   r::AbstractVector{Float64},
                   n::AbstractVector{Float64},
                   X::AbstractMatrix{Float64},
                   beta::AbstractVector{Float64},
                   M::Float64)::Float64
   m = length(k)
   p = length(beta)
   ll = 0.0
   lgM = _loggamma(M)
   @inbounds for j in 1:m
      eta = 0.0
      for c in 1:p
         eta += X[j, c] * beta[c]
      end
      mu = clamp(_logistic(eta), _MU_EPS, 1.0 - _MU_EPS)
      a = M * mu
      b = M * (1.0 - mu)
      ll += _loggamma(k[j] + a) + _loggamma(r[j] + b) - _loggamma(n[j] + M) +
            lgM - _loggamma(a) - _loggamma(b)
   end
   return ll
end

# Fit beta for fixed precision M using damped Fisher scoring (IRLS) with tiny ridge
# regularization for numerical stability under complete separation (e.g. 0 vs 1 PSI).
# Returns (beta, unpenalized_ll, fisher_info, converged).
function fit_bb_beta(k::AbstractVector{Float64},
                     r::AbstractVector{Float64},
                     n::AbstractVector{Float64},
                     X::AbstractMatrix{Float64},
                     M::Float64;
                     beta_init::Union{Nothing, Vector{Float64}}=nothing,
                     ridge::Float64=1e-4,
                     maxiter::Int=40,
                     tol::Float64=1e-8)
   m, p = size(X)
   beta = beta_init === nothing ? zeros(Float64, p) : copy(beta_init)
   if beta_init === nothing && p >= 1
      # Initialize intercept at depth-weighted mean PSI
      sk = sum(k)
      sn = sum(n)
      if sn > 0.0
         beta[1] = _logit(clamp(sk / sn, 1e-4, 1.0 - 1e-4))
      end
   end

   g = Vector{Float64}(undef, p)
   H = Matrix{Float64}(undef, p, p)
   beta_try = Vector{Float64}(undef, p)

   @inline function pen_val(b::Vector{Float64})::Float64
      pv = 0.0
      @inbounds for c in 1:p
         lam = c == 1 ? (ridge * 0.01) : ridge
         pv += 0.5 * lam * b[c] * b[c]
      end
      return pv
   end

   ll_cur = bb_loglik(k, r, n, X, beta, M)
   obj_cur = ll_cur - pen_val(beta)
   converged = false

   for iter in 1:maxiter
      fill!(g, 0.0)
      fill!(H, 0.0)
      @inbounds for j in 1:m
         eta = 0.0
         for c in 1:p
            eta += X[j, c] * beta[c]
         end
         mu = clamp(_logistic(eta), _MU_EPS, 1.0 - _MU_EPS)
         mu1 = mu * (1.0 - mu)
         a = M * mu
         b = M * (1.0 - mu)

         # Exact score d(ll_j)/d(eta_j)
         score_j = M * mu1 * (_digamma(k[j] + a) - _digamma(a) -
                              _digamma(r[j] + b) + _digamma(b))
         # Positive-definite Fisher scoring weight
         t_term = (_trigamma(a) - _trigamma(k[j] + a)) +
                  (_trigamma(b) - _trigamma(r[j] + b))
         w_j = max(M * M * mu1 * mu1 * t_term, 1e-12)

         for c in 1:p
            xjc = X[j, c]
            g[c] += xjc * score_j
            wx = w_j * xjc
            for d in c:p
               H[c, d] += wx * X[j, d]
            end
         end
      end

      # Mirror upper triangle & apply ridge penalty
      @inbounds for c in 1:p
         lam = c == 1 ? (ridge * 0.01) : ridge
         g[c] -= lam * beta[c]
         H[c, c] += lam
         for d in (c + 1):p
            H[d, c] = H[c, d]
         end
      end

      # Solve H * step = g with Cholesky fallback to damped diagonal
      step = Vector{Float64}(undef, p)
      ch = cholesky(Symmetric(H); check=false)
      if issuccess(ch)
         step .= ch \ g
      else
         H_reg = copy(H)
         diag_bump = 1e-6
         solved = false
         for _ in 1:8
            for c in 1:p
               H_reg[c, c] = H[c, c] + diag_bump
            end
            ch2 = cholesky(Symmetric(H_reg); check=false)
            if issuccess(ch2)
               step .= ch2 \ g
               solved = true
               break
            end
            diag_bump *= 10.0
         end
         if !solved
            break
         end
      end

      # Clamp extreme Newton step for stability
      max_step = maximum(abs, step)
      if max_step > 5.0
         scale = 5.0 / max_step
         @inbounds for c in 1:p
            step[c] *= scale
         end
      end

      # Backtracking line search ensuring monotone ascent of penalized log-likelihood
      alpha = 1.0
      accepted = false
      ll_new = ll_cur
      obj_new = obj_cur
      for _ in 1:20
         @inbounds for c in 1:p
            beta_try[c] = beta[c] + alpha * step[c]
         end
         ll_new = bb_loglik(k, r, n, X, beta_try, M)
         obj_new = ll_new - pen_val(beta_try)
         if isfinite(obj_new) && obj_new >= obj_cur - 1e-12
            accepted = true
            break
         end
         alpha *= 0.5
      end

      if !accepted
         converged = true
         break
      end

      max_delta = 0.0
      @inbounds for c in 1:p
         d = abs(beta_try[c] - beta[c])
         d > max_delta && (max_delta = d)
         beta[c] = beta_try[c]
      end
      d_obj = abs(obj_new - obj_cur)
      ll_cur = ll_new
      obj_cur = obj_new

      if max_delta < tol || d_obj < tol * (1.0 + abs(obj_cur))
         converged = true
         break
      end
   end

   # Recompute Fisher info matrix H at converged beta (without ridge for Cox-Reid logdet)
   fill!(H, 0.0)
   @inbounds for j in 1:m
      eta = 0.0
      for c in 1:p
         eta += X[j, c] * beta[c]
      end
      mu = clamp(_logistic(eta), _MU_EPS, 1.0 - _MU_EPS)
      mu1 = mu * (1.0 - mu)
      a = M * mu
      b = M * (1.0 - mu)
      t_term = (_trigamma(a) - _trigamma(k[j] + a)) +
               (_trigamma(b) - _trigamma(r[j] + b))
      w_j = max(M * M * mu1 * mu1 * t_term, 1e-12)
      for c in 1:p
         wx = w_j * X[j, c]
         for d in c:p
            H[c, d] += wx * X[j, d]
         end
      end
   end
   @inbounds for c in 1:p
      for d in (c + 1):p
         H[d, c] = H[c, d]
      end
   end

   return beta, ll_cur, H, converged
end

# Cox-Reid Adjusted Profile Log-Likelihood at theta = log(M):
#   APL(theta) = ll(beta_hat(theta), exp(theta)) - 0.5 * logdet(X' W X)
function bb_cox_reid_apl(k::AbstractVector{Float64},
                         r::AbstractVector{Float64},
                         n::AbstractVector{Float64},
                         X::AbstractMatrix{Float64},
                         theta::Float64;
                         beta_init::Union{Nothing, Vector{Float64}}=nothing,
                         ridge::Float64=1e-4)
   M = exp(theta)
   beta, ll, H, _ = fit_bb_beta(k, r, n, X, M; beta_init=beta_init, ridge=ridge)
   p = size(X, 2)
   # Small jitter on H diagonal only if near-singular (complete separation)
   H_cr = copy(H)
   @inbounds for c in 1:p
      H_cr[c, c] += 1e-8
   end
   ch = cholesky(Symmetric(H_cr); check=false)
   cr_pen = issuccess(ch) ? 0.5 * logdet(ch) : 0.0
   return ll - cr_pen, beta
end

# 1D maximization on [lo, hi] via coarse grid + Brent/golden-section refinement.
# Warm-starts beta across evaluations for speed.
function _maximize_theta(f, lo::Float64, hi::Float64; ngrid::Int=10, tol::Float64=1e-3)
   # Coarse grid to guarantee finding global basin (APL can be flat at extremes)
   best_t = lo
   best_v = -Inf
   best_idx = 1
   grid = collect(range(lo, hi, length=ngrid))
   vals = Vector{Float64}(undef, ngrid)
   for i in 1:ngrid
      v = f(grid[i])
      vals[i] = v
      if v > best_v
         best_v = v
         best_t = grid[i]
         best_idx = i
      end
   end

   a = best_idx > 1 ? grid[best_idx - 1] : lo
   b = best_idx < ngrid ? grid[best_idx + 1] : hi

   # Golden-section refinement on [a, b]
   invphi = 0.6180339887498949
   invphi2 = 0.3819660112501051
   h = b - a
   c = a + invphi2 * h
   d = a + invphi * h
   fc = f(c)
   fd = f(d)
   for _ in 1:35
      if (b - a) < tol
         break
      end
      if fc > fd
         b = d
         d = c
         fd = fc
         h = b - a
         c = a + invphi2 * h
         fc = f(c)
      else
         a = c
         c = d
         fc = fd
         h = b - a
         d = a + invphi * h
         fd = f(d)
      end
   end

   t_opt = 0.5 * (a + b)
   v_opt = f(t_opt)
   if best_v > v_opt
      t_opt = best_t
      v_opt = best_v
   end
   return t_opt, v_opt
end

# Estimate raw dispersion theta_hat = log(M_hat) and its local curvature (Fisher info for theta)
function estimate_raw_theta(k::AbstractVector{Float64},
                            r::AbstractVector{Float64},
                            n::AbstractVector{Float64},
                            X::AbstractMatrix{Float64};
                            ridge::Float64=1e-4)
   beta_warm = Ref{Union{Nothing, Vector{Float64}}}(nothing)
   function eval_apl(th::Float64)::Float64
      val, b = bb_cox_reid_apl(k, r, n, X, th; beta_init=beta_warm[], ridge=ridge)
      beta_warm[] = b
      return val
   end

   th_hat, apl_hat = _maximize_theta(eval_apl, _THETA_MIN, _THETA_MAX; ngrid=9, tol=2e-3)

   # Estimate sampling variance 1 / (-APL''(th_hat)) via 3-point finite difference if interior
   is_interior = (_THETA_MIN + 0.05) < th_hat < (_THETA_MAX - 0.05)
   samp_var = NaN
   if is_interior
      h = 0.1
      v_lo = eval_apl(th_hat - h)
      v_hi = eval_apl(th_hat + h)
      curv = -(v_hi - 2.0 * apl_hat + v_lo) / (h * h)
      if curv > 1e-4
         samp_var = 1.0 / curv
      end
   end
   return th_hat, is_interior, samp_var
end

# Estimate Empirical Bayes MAP shrunk theta = log(M) under Normal(prior_mu, prior_var)
function estimate_eb_theta(k::AbstractVector{Float64},
                           r::AbstractVector{Float64},
                           n::AbstractVector{Float64},
                           X::AbstractMatrix{Float64},
                           prior_mu::Float64,
                           prior_var::Float64;
                           ridge::Float64=1e-4)
   beta_warm = Ref{Union{Nothing, Vector{Float64}}}(nothing)
   inv_2v = 0.5 / prior_var
   function eval_post(th::Float64)::Float64
      val, b = bb_cox_reid_apl(k, r, n, X, th; beta_init=beta_warm[], ridge=ridge)
      beta_warm[] = b
      diff = th - prior_mu
      return val - diff * diff * inv_2v
   end

   th_map, _ = _maximize_theta(eval_post, _THETA_MIN, _THETA_MAX; ngrid=9, tol=2e-3)
   return th_map
end

# Estimate genome-wide Empirical Bayes prior Normal(mu_theta, sigma_theta^2) across raw event estimates.
# Applies the exact log-chi2 median bias correction log(median(Chisq(d_res)) / d_res) so the
# prior mean is unbiased for small residual degrees of freedom (Smyth 2004 / edgeR).
function fit_eb_dispersion_prior(raw_thetas::AbstractVector{Float64},
                                 is_interior::AbstractVector{Bool},
                                 samp_vars::AbstractVector{Float64},
                                 res_dfs::AbstractVector{Int}=Int[])
   interior_th = Float64[]
   interior_sv = Float64[]
   interior_df = Float64[]
   all_valid   = Float64[]
   has_dfs = length(res_dfs) == length(raw_thetas)
   @inbounds for i in 1:length(raw_thetas)
      t = raw_thetas[i]
      isfinite(t) || continue
      push!(all_valid, t)
      if is_interior[i]
         push!(interior_th, t)
         if isfinite(samp_vars[i])
            push!(interior_sv, samp_vars[i])
         end
         if has_dfs && res_dfs[i] >= 1
            push!(interior_df, Float64(res_dfs[i]))
         end
      end
   end

   # Finite-df log-chi2 median bias correction: E/median of -log(chi2_d / d) > 0 for small d
   bias_corr = 0.0
   if !isempty(interior_df)
      med_df = max(1.0, Distributions.Statistics.median(interior_df))
      bias_corr = log(quantile(Chisq(med_df), 0.5) / med_df)
   end

   if length(interior_th) >= 20
      raw_med  = Distributions.Statistics.median(interior_th)
      mu_prior = clamp(raw_med + bias_corr, _THETA_MIN, _THETA_MAX)
      mad_val  = Distributions.Statistics.median(abs.(interior_th .- raw_med))
      obs_sd   = 1.4826 * mad_val
      med_sv   = isempty(interior_sv) ? 0.0 : Distributions.Statistics.median(interior_sv)
      var_prior = clamp(obs_sd * obs_sd - med_sv, 0.08, 9.0)
      return mu_prior, var_prior, length(interior_th)
   elseif length(all_valid) >= 10
      raw_med  = Distributions.Statistics.median(all_valid)
      mu_prior = clamp(raw_med + bias_corr, _THETA_MIN, _THETA_MAX)
      mad_val  = Distributions.Statistics.median(abs.(all_valid .- raw_med))
      obs_sd   = 1.4826 * mad_val
      var_prior = clamp(obs_sd * obs_sd, 0.25, 9.0)
      return mu_prior, var_prior, length(all_valid)
   else
      # Weakly informative fallback prior (M ≈ 20 => phi ≈ 0.047)
      return 3.0, 4.0, length(all_valid)
   end
end

# Benjamini-Hochberg FDR adjustment (handles NaN cleanly by preserving NaN)
function bh_fdr(pvals::AbstractVector{Float64})::Vector{Float64}
   n = length(pvals)
   fdr = fill(NaN, n)
   valid_idx = Int[]
   sizehint!(valid_idx, n)
   @inbounds for i in 1:n
      p = pvals[i]
      if isfinite(p)
         push!(valid_idx, i)
      end
   end
   m = length(valid_idx)
   m == 0 && return fdr

   sort!(valid_idx, by = i -> pvals[i])
   running_min = 1.0
   @inbounds for rank in m:-1:1
      idx = valid_idx[rank]
      q = pvals[idx] * (Float64(m) / Float64(rank))
      if q < running_min
         running_min = q
      end
      fdr[idx] = clamp(running_min, 0.0, 1.0)
   end
   return fdr
end

# ============================================================================
# Design Matrix & Formula Parsing (DESeq2-style: "~ batch + donor + condition")
# ============================================================================

struct DesignSpec
   files::Vector{String}
   sample_names::Vector{String}
   X::Matrix{Float64}
   col_names::Vector{String}
   test_cols::Vector{Int}           # column indices in X being tested (dropped in null model)
   is_two_group::Bool               # true when comparing two explicit levels (A vs B)
   group_a_mask::BitVector          # samples in group A (numerator / first group)
   group_b_mask::BitVector          # samples in group B (denominator / reference group)
   group_a_label::String
   group_b_label::String
   test_label::String
   is_continuous_test::Bool         # true when testing a single numeric covariate
   delta_sign::Float64              # +1.0 for Psi_A - Psi_B
end

# Build simple 2-group design matrix for classic `-a` / `-b` mode.
# Convention matches original whippet-delta.jl: DeltaPsi = Psi_A - Psi_B.
# Design: column 1 = intercept (group B baseline), column 2 = indicator(sample in A).
# Thus beta[2] > 0 means Psi_A > Psi_B (DeltaPsi > 0)!
function build_two_group_design(files_a::Vector{String},
                                files_b::Vector{String};
                                label_a::String="A",
                                label_b::String="B")::DesignSpec
   na = length(files_a)
   nb = length(files_b)
   ns = na + nb
   files = vcat(files_a, files_b)
   snames = [basename(f) for f in files]
   X = zeros(Float64, ns, 2)
   @inbounds for i in 1:ns
      X[i, 1] = 1.0
      X[i, 2] = i <= na ? 1.0 : 0.0
   end
   mask_a = BitVector([i <= na for i in 1:ns])
   mask_b = BitVector([i > na for i in 1:ns])
   return DesignSpec(files, snames, X, ["(Intercept)", "$(label_a)_vs_$(label_b)"],
                     [2], true, mask_a, mask_b, label_a, label_b,
                     "$(label_a)_vs_$(label_b)", false, 1.0)
end

# Parse a formula string like "~ batch + donor + condition" or "~ batch * condition"
# into an ordered list of additive terms (each term is a Vector{String} of interacting variable names).
function parse_formula_terms(formula::AbstractString)::Vector{Vector{String}}
   s = strip(String(formula))
   if startswith(s, "~")
      s = strip(s[2:end])
   end
   isempty(s) && error("Empty formula: '$formula'")

   raw_parts = [String(strip(p)) for p in split(s, '+') if !isempty(strip(p))]
   terms = Vector{Vector{String}}()
   seen = Set{Vector{String}}()

   function add_term!(vars::Vector{String})
      canon = sort(unique(vars))
      if !(canon in seen)
         push!(seen, canon)
         push!(terms, vars)
      end
   end

   for part in raw_parts
      part == "1" && continue
      if occursin('*', part)
         factors = [String(strip(x)) for x in split(part, '*') if !isempty(strip(x))]
         # Expand a*b -> a, b, a:b (supports 2-way or k-way interactions)
         k = length(factors)
         for mask in 1:((1 << k) - 1)
            sub = String[]
            for bit in 0:(k - 1)
               if ((mask >> bit) & 1) == 1
                  push!(sub, factors[bit + 1])
               end
            end
            add_term!(sub)
         end
      elseif occursin(':', part)
         factors = [String(strip(x)) for x in split(part, ':') if !isempty(strip(x))]
         add_term!(factors)
      else
         add_term!(String[part])
      end
   end
   isempty(terms) && error("Formula '$formula' contains no predictor terms!")
   return terms
end

# Try to treat a DataFrame column as numeric continuous covariate; returns (is_num, Vector{Float64})
function _try_numeric_col(col)
   n = length(col)
   vals = Vector{Float64}(undef, n)
   for i in 1:n
      v = col[i]
      if ismissing(v)
         return false, vals
      elseif v isa Number && !(v isa Bool)
         vals[i] = Float64(v)
      else
         p = tryparse(Float64, strip(string(v)))
         p === nothing && return false, vals
         vals[i] = p
      end
   end
   return true, vals
end

# Encode a single variable column into design columns:
# - Numeric (non-categorical) -> 1 column (raw float values)
# - Categorical -> (K-1) dummy indicator columns with reference level ref_level
function _encode_variable(df::DataFrame, varname::String; ref_level::Union{Nothing, String}=nothing,
                          force_categorical::Bool=false)
   hasproperty(df, Symbol(varname)) || error("Variable '$varname' not found in design table columns: $(names(df))")
   col = df[!, Symbol(varname)]
   n = nrow(df)

   if !force_categorical
      is_num, num_vals = _try_numeric_col(col)
      if is_num
         return reshape(num_vals, n, 1), [varname], false, String[]
      end
   end

   str_vals = [strip(string(v)) for v in col]
   levels = sort(unique(str_vals))
   length(levels) >= 2 || error("Categorical variable '$varname' has fewer than 2 levels: $levels")

   if ref_level !== nothing
      ref_level in levels || error("Reference level '$ref_level' not found in '$varname' (available: $levels)")
      other_levels = [lv for lv in levels if lv != ref_level]
      ordered_levels = vcat([ref_level], other_levels)
   else
      ordered_levels = levels
      other_levels = levels[2:end]
   end

   M = zeros(Float64, n, length(other_levels))
   cnames = String[]
   for (j, lv) in enumerate(other_levels)
      push!(cnames, "$(varname):$(lv)")
      for i in 1:n
         if str_vals[i] == lv
            M[i, j] = 1.0
         end
      end
   end
   return M, cnames, true, ordered_levels
end

# Load a CSV/TSV design table and build the full DesignSpec from formula & test specification.
# Expected columns: `file` (path to .psi/.psi.gz, resolved relative to `dir` if not absolute),
# optional `sample`, plus any predictor/covariate columns referenced in `formula`.
# `test_spec` can be:
#   - "" (default: tests the last term in `formula`)
#   - "condition" (tests all columns belonging to `condition`)
#   - "condition:treated_vs_ctrl" (sets `ctrl` as reference level and tests `treated` vs `ctrl`,
#     reporting Psi_A = treated, Psi_B = ctrl, DeltaPsi = Psi(treated) - Psi(ctrl))
function build_design_from_table(design_path::String,
                                 formula::String;
                                 test_spec::String="",
                                 dir::String=".",
                                 file_col::String="file")::DesignSpec
   df = DataFrame(CSV.File(design_path))
   fcol_sym = Symbol(file_col)
   hasproperty(df, fcol_sym) || error("Design table '$design_path' must have a '$file_col' column (use --file-col to specify another column name; available: $(names(df))).")
   ns = nrow(df)
   ns >= 2 || error("Design table must have at least 2 samples.")

   files = String[]
   for f in df[!, fcol_sym]
      fs = strip(string(f))
      full = isabspath(fs) ? fs : joinpath(dir, fs)
      if !isfile(full) && isfile(joinpath(dirname(abspath(design_path)), fs))
         full = joinpath(dirname(abspath(design_path)), fs)
      end
      isfile(full) || error("PSI file not found: '$full' (from design table entry '$fs')")
      push!(files, full)
   end
   snames = hasproperty(df, :sample) ? [strip(string(s)) for s in df[!, :sample]] : [basename(f) for f in files]

   terms = parse_formula_terms(formula)

   # Parse test_spec: either "var", "varA:varB", or "var:levelA_vs_levelB"
   test_var = ""
   contrast_num = nothing
   contrast_den = nothing
   if isempty(strip(test_spec))
      test_var = join(terms[end], ":")
   elseif occursin("_vs_", test_spec) && occursin(':', test_spec)
      idx = findfirst(isequal(':'), test_spec)
      test_var = String(strip(test_spec[1:(idx - 1)]))
      rhs = String(strip(test_spec[(idx + 1):end]))
      parts = split(rhs, "_vs_")
      length(parts) == 2 || error("Invalid contrast specification '$test_spec'. Expected 'column:LevelA_vs_LevelB'.")
      contrast_num = String(strip(parts[1]))
      contrast_den = String(strip(parts[2]))
   else
      test_var = String(strip(test_spec))
   end

   # Cache per-variable encodings
   var_mats = Dict{String, Matrix{Float64}}()
   var_cnames = Dict{String, Vector{String}}()
   var_iscat = Dict{String, Bool}()
   var_levels = Dict{String, Vector{String}}()

   all_vars = unique(vcat(terms...))
   for v in all_vars
      ref_lv = (v == test_var && contrast_den !== nothing) ? contrast_den : nothing
      force_cat = (v == test_var && contrast_den !== nothing)
      M, cn, iscat, lvs = _encode_variable(df, v; ref_level=ref_lv, force_categorical=force_cat)
      var_mats[v] = M
      var_cnames[v] = cn
      var_iscat[v] = iscat
      var_levels[v] = lvs
   end

   # Assemble full design matrix X with intercept
   cols_list = Matrix{Float64}[ones(Float64, ns, 1)]
   col_names = String["(Intercept)"]
   term_col_map = Dict{String, Vector{Int}}()
   cur_col = 1

   for term in terms
      tname = join(term, ":")
      if length(term) == 1
         v = term[1]
         M = var_mats[v]
         cns = var_cnames[v]
      else
         # Elementwise column products across all interacting variables
         M = var_mats[term[1]]
         cns = var_cnames[term[1]]
         for k in 2:length(term)
            M2 = var_mats[term[k]]
            cns2 = var_cnames[term[k]]
            M_new = Matrix{Float64}(undef, ns, size(M, 2) * size(M2, 2))
            cns_new = String[]
            col_i = 0
            for c1 in 1:size(M, 2)
               for c2 in 1:size(M2, 2)
                  col_i += 1
                  @inbounds for row in 1:ns
                     M_new[row, col_i] = M[row, c1] * M2[row, c2]
                  end
                  push!(cns_new, "$(cns[c1]):$(cns2[c2])")
               end
            end
            M = M_new
            cns = cns_new
         end
      end
      ncols = size(M, 2)
      push!(cols_list, M)
      append!(col_names, cns)
      term_col_map[tname] = collect((cur_col + 1):(cur_col + ncols))
      cur_col += ncols
   end

   X = hcat(cols_list...)

   # Identify tested columns and whether this is a 2-group contrast
   haskey(term_col_map, test_var) || error("Test term '$test_var' not found in formula terms: $(collect(keys(term_col_map)))")
   test_cols = term_col_map[test_var]

   is_two_group = false
   mask_a = falses(ns)
   mask_b = falses(ns)
   label_a = "A"
   label_b = "B"
   is_cont = false
   tlabel = test_var

   if haskey(var_iscat, test_var)
      if !var_iscat[test_var]
         is_cont = true
      else
         str_vals = [strip(string(v)) for v in df[!, Symbol(test_var)]]
         lvs = var_levels[test_var]
         if contrast_num !== nothing && contrast_den !== nothing
            contrast_num in lvs || error("Contrast level '$contrast_num' not found in '$test_var' (available: $lvs)")
            # Find the specific column for contrast_num
            target_col_name = "$(test_var):$(contrast_num)"
            cidx = findfirst(==(target_col_name), col_names)
            cidx === nothing && error("Could not locate contrast column '$target_col_name'")
            test_cols = [cidx]
            is_two_group = true
            mask_a = BitVector([s == contrast_num for s in str_vals])
            mask_b = BitVector([s == contrast_den for s in str_vals])
            label_a = contrast_num
            label_b = contrast_den
            tlabel = "$(test_var):$(contrast_num)_vs_$(contrast_den)"
         elseif length(lvs) == 2
            # Binary categorical variable: lvs[1] = reference (B), lvs[2] = non-reference (A)
            is_two_group = true
            mask_a = BitVector([s == lvs[2] for s in str_vals])
            mask_b = BitVector([s == lvs[1] for s in str_vals])
            label_a = lvs[2]
            label_b = lvs[1]
            tlabel = "$(test_var):$(lvs[2])_vs_$(lvs[1])"
         end
      end
   end

   return DesignSpec(files, snames, X, col_names, test_cols, is_two_group,
                     mask_a, mask_b, label_a, label_b, tlabel, is_cont, 1.0)
end

# ============================================================================
# Event Streaming, Per-Event Collinearity Handling, Parallel GLM & Output
# ============================================================================

struct GLMEventData
   gene::String
   node::String
   coord::String
   strand::String
   etype::String
   complex::Int
   entropy::Float64
   valid_idx::Vector{Int}
   psi::Vector{Float64}
   reads::Vector{Float64}
   samp_a::Int
   samp_b::Int
   mean_reads_a::Float64
   mean_reads_b::Float64
end

struct GLMResult
   tested::Bool
   psi_a::Float64
   psi_b::Float64
   delta_psi::Float64
   log_odds::Float64
   se::Float64
   dispersion::Float64
   lrt::Float64
   df::Int
   pvalue::Float64
end

const _EMPTY_GLM_RESULT = GLMResult(false, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 0, NaN)

# Build full-rank reduced/full sub-matrices for an event after subsetting to `valid_idx`.
# Always retains column 1 ((Intercept)) first, drops collinear nuisance columns
# (e.g. a batch level absent in `valid_idx`), and returns `nothing` if the tested
# contrast itself becomes unidentifiable.
function _build_event_matrices(X::Matrix{Float64},
                               valid_idx::Vector{Int},
                               test_cols::Vector{Int})
   m = length(valid_idx)
   p_all = size(X, 2)
   test_set = Set{Int}(test_cols)
   nuis_cols = [c for c in 1:p_all if !(c in test_set)]

   function _adds_rank(base::Matrix{Float64}, col::Vector{Float64})::Bool
      cand = size(base, 2) == 0 ? reshape(col, m, 1) : hcat(base, col)
      F = qr(cand, ColumnNorm())
      R = F.R
      rmax = maximum(abs, diag(R))
      tol = max(m, size(cand, 2)) * eps(Float64) * max(rmax, 1.0)
      rnk = 0
      for i in 1:min(size(R)...)
         if abs(R[i, i]) > tol
            rnk += 1
         end
      end
      return rnk > size(base, 2)
   end

   # 1. Extract linearly independent nuisance columns in deterministic order (intercept first)
   X0 = zeros(Float64, m, 0)
   for nc in nuis_cols
      col = X[valid_idx, nc]
      if _adds_rank(X0, col)
         X0 = size(X0, 2) == 0 ? reshape(col, m, 1) : hcat(X0, col)
      end
   end
   p0 = size(X0, 2)

   # 2. Keep tested columns that add independent rank beyond X0
   keep_test = Int[]
   cur_X = X0
   for tc in test_cols
      col = X[valid_idx, tc]
      if _adds_rank(cur_X, col)
         push!(keep_test, tc)
         cur_X = size(cur_X, 2) == 0 ? reshape(col, m, 1) : hcat(cur_X, col)
      end
   end

   # Must have at least 1 identifiable tested df and enough samples (m >= size(cur_X, 2))
   if isempty(keep_test) || m < size(cur_X, 2)
      return nothing
   end

   # Tested columns sit at the end of cur_X: indices (p0 + 1) : (p0 + length(keep_test))
   return cur_X, X0, p0, keep_test
end

# Read all `.psi.gz` files in lockstep, keeping canonical events that pass read/sample filters.
function read_psi_cohort(design::DesignSpec;
                         min_reads::Int=5,
                         min_samp::Int=1,
                         whitelist_events::Set{String}=_GLM_WHITELIST_EVENTS)::Vector{GLMEventData}
   ns = length(design.files)
   streams = open_streams(design.files)
   events = GLMEventData[]
   min_r = Float64(min_reads)

   # Pre-allocated scratch per row
   row_psi   = Vector{Float64}(undef, ns)
   row_reads = Vector{Float64}(undef, ns)
   row_valid = falses(ns)

   try
      while true
         ref_gene   = ""
         ref_node   = ""
         ref_coord  = ""
         ref_strand = ""
         ref_type   = ""
         hit_eof    = false
         max_comp   = 0
         max_ent    = 0.0

         for s in 1:ns
            row_valid[s] = false
            row_psi[s] = NaN
            row_reads[s] = 0.0

            # Advance stream s until a whitelisted canonical event line or EOF
            while true
               if eof(streams[s])
                  hit_eof = true
                  break
               end
               line = readline(streams[s])
               isempty(line) && continue
               fields = split(line, '\t')
               length(fields) >= 11 || continue
               etype = String(fields[5])
               etype in whitelist_events || continue

               if s == 1
                  ref_gene   = String(fields[1])
                  ref_node   = String(fields[2])
                  ref_coord  = String(fields[3])
                  ref_strand = String(fields[4])
                  ref_type   = etype
               else
                  if fields[1] != ref_gene || fields[2] != ref_node || fields[5] != ref_type
                     error("Mismatched events across PSI files! Expected ($ref_gene, $ref_node, $ref_type) " *
                           "from $(basename(design.files[1])), got ($(fields[1]), $(fields[2]), $(fields[5])) " *
                           "in $(basename(design.files[s])). Ensure all samples were quantified with the same index.")
                  end
               end

               psi_str = fields[6]
               num_str = fields[9]
               if psi_str != "NA" && num_str != "NA"
                  psi_v = tryparse(Float64, psi_str)
                  num_v = tryparse(Float64, num_str)
                  if psi_v !== nothing && num_v !== nothing && psi_v >= 0.0 && num_v >= min_r
                     row_psi[s] = clamp(psi_v, 0.0, 1.0)
                     row_reads[s] = num_v
                     row_valid[s] = true
                     comp_str = fields[10]
                     if comp_str != "NA" && startswith(comp_str, COMPLEX_CHAR)
                        cv = tryparse(Int, comp_str[2:end])
                        cv !== nothing && cv > max_comp && (max_comp = cv)
                     end
                     ent_str = fields[11]
                     if ent_str != "NA"
                        ev = tryparse(Float64, ent_str)
                        ev !== nothing && ev > max_ent && (max_ent = ev)
                     end
                  end
               end
               break
            end
            hit_eof && break
         end
         hit_eof && break

         # Apply sample-count filter
         if design.is_two_group
            na = 0; nb = 0
            sum_ra = 0.0; sum_rb = 0.0
            @inbounds for s in 1:ns
               row_valid[s] || continue
               if design.group_a_mask[s]
                  na += 1
                  sum_ra += row_reads[s]
               elseif design.group_b_mask[s]
                  nb += 1
                  sum_rb += row_reads[s]
               end
            end
            if na >= min_samp && nb >= min_samp
               vidx = Int[]
               vpsi = Float64[]
               vrds = Float64[]
               @inbounds for s in 1:ns
                  if row_valid[s]
                     push!(vidx, s)
                     push!(vpsi, row_psi[s])
                     push!(vrds, row_reads[s])
                  end
               end
               push!(events, GLMEventData(ref_gene, ref_node, ref_coord, ref_strand, ref_type,
                                          max_comp, max_ent, vidx, vpsi, vrds,
                                          na, nb, sum_ra / na, sum_rb / nb))
            end
         else
            nv = 0
            sum_r = 0.0
            @inbounds for s in 1:ns
               if row_valid[s]
                  nv += 1
                  sum_r += row_reads[s]
               end
            end
            if nv >= max(min_samp, 2)
               vidx = Int[]
               vpsi = Float64[]
               vrds = Float64[]
               @inbounds for s in 1:ns
                  if row_valid[s]
                     push!(vidx, s)
                     push!(vpsi, row_psi[s])
                     push!(vrds, row_reads[s])
                  end
               end
               push!(events, GLMEventData(ref_gene, ref_node, ref_coord, ref_strand, ref_type,
                                          max_comp, max_ent, vidx, vpsi, vrds,
                                          nv, 0, sum_r / nv, 0.0))
            end
         end
      end
   finally
      for st in streams
         try close(st) catch end
      end
   end

   return events
end

# Fit full + null Beta-Binomial GLM and compute LRT + marginal effect sizes for a single event.
# When `df_denom` is finite (estimated dispersion with residual + EB prior degrees of freedom),
# uses the finite-sample F-distribution F(df, df_denom) on LRT/df (matching edgeR glmQLFTest /
# limma empirical Bayes F-test); when `df_denom == Inf` (fixed known dispersion), reduces to Chisq(df).
function _test_single_event(ev::GLMEventData,
                            design::DesignSpec,
                            theta_use::Float64;
                            df_denom::Float64=Inf,
                            ridge::Float64=1e-4)::GLMResult
   mats = _build_event_matrices(design.X, ev.valid_idx, design.test_cols)
   mats === nothing && return _EMPTY_GLM_RESULT
   X_full, X_null, p0, keep_test = mats
   p_full = size(X_full, 2)
   df = p_full - p0
   df >= 1 || return _EMPTY_GLM_RESULT

   m = length(ev.valid_idx)
   k = Vector{Float64}(undef, m)
   r = Vector{Float64}(undef, m)
   n = ev.reads
   @inbounds for j in 1:m
      kj = ev.psi[j] * n[j]
      k[j] = kj
      r[j] = n[j] - kj
   end

   M = exp(theta_use)
   phi = 1.0 / (1.0 + M)

   # Fit full model at shrunk M
   beta_full, ll_full, H_full, _ = fit_bb_beta(k, r, n, X_full, M; ridge=ridge)

   # Fit null model at the same M
   if p0 > 0
      beta_null_init = beta_full[1:p0]
      _, ll_null, _, _ = fit_bb_beta(k, r, n, X_null, M; beta_init=beta_null_init, ridge=ridge)
   else
      # Edge case: no intercept/nuisance columns (null mu = 0.5)
      ll_null = bb_loglik(k, r, n, zeros(Float64, m, 1), [0.0], M)
   end

   lrt = max(0.0, 2.0 * (ll_full - ll_null))
   if !isfinite(lrt)
      pval = NaN
   elseif isfinite(df_denom) && df_denom > 0.0
      pval = ccdf(FDist(Float64(df), df_denom), lrt / Float64(df))
   else
      pval = ccdf(Chisq(df), lrt)
   end

   # Extract coefficient & SE for 1-df tests (first tested column)
   log_odds = NaN
   se = NaN
   if df == 1
      t_col_pos = p_full
      log_odds = beta_full[t_col_pos]
      H_inv_mat = copy(H_full)
      @inbounds for c in 1:p_full
         H_inv_mat[c, c] += (c == 1 ? ridge * 0.01 : ridge)
      end
      ch = cholesky(Symmetric(H_inv_mat); check=false)
      if issuccess(ch)
         e_vec = zeros(Float64, p_full)
         e_vec[t_col_pos] = 1.0
         v_col = ch \ e_vec
         v_diag = v_col[t_col_pos]
         if v_diag >= 0.0 && isfinite(v_diag)
            se = sqrt(v_diag)
         end
      end
   end

   # Compute covariate-adjusted marginal Psi_A, Psi_B, and DeltaPsi
   psi_a = NaN
   psi_b = NaN
   delta_psi = NaN

   if design.is_two_group && df == 1
      t_col_pos = p_full
      # Marginal standardization across valid samples belonging to Group A or Group B:
      # predict with tested contrast set to 1 (Group A) vs 0 (Group B), holding covariates fixed.
      # Restricting to A/B samples ensures other factor levels (in 3+ level contrasts) do not leak into eta_base.
      sum_mu_a = 0.0
      sum_mu_b = 0.0
      n_ab = 0
      @inbounds for j in 1:m
         s_orig = ev.valid_idx[j]
         (design.group_a_mask[s_orig] || design.group_b_mask[s_orig]) || continue
         eta_base = 0.0
         for c in 1:p0
            eta_base += X_full[j, c] * beta_full[c]
         end
         sum_mu_b += _logistic(eta_base)
         sum_mu_a += _logistic(eta_base + beta_full[t_col_pos])
         n_ab += 1
      end
      if n_ab > 0
         psi_a = sum_mu_a / n_ab
         psi_b = sum_mu_b / n_ab
         delta_psi = psi_a - psi_b
      end
   else
      # Continuous covariate or multi-level factor: report depth-weighted mean PSI
      sk = sum(k)
      sn = sum(n)
      psi_a = sn > 0.0 ? (sk / sn) : NaN
      if design.is_continuous_test && df == 1 && isfinite(log_odds)
         # Average marginal effect d(Psi)/dx = mean_j(mu_j * (1 - mu_j) * beta)
         sum_slope = 0.0
         @inbounds for j in 1:m
            eta = 0.0
            for c in 1:p_full
               eta += X_full[j, c] * beta_full[c]
            end
            mu = _logistic(eta)
            sum_slope += mu * (1.0 - mu) * log_odds
         end
         delta_psi = sum_slope / m
      end
   end

   return GLMResult(true, psi_a, psi_b, delta_psi, log_odds, se, phi, lrt, df, pval)
end

# Returns true iff an event is polymorphic (non-constitutive) and carries statistical
# information about Beta-Binomial overdispersion phi. Constitutive/boundary events
# (Psi == 0.0 or 1.0 across all samples) have a flat Beta-Binomial log-likelihood in M,
# causing Cox-Reid -0.5 * logdet(X'WX) to create an artificial interior peak near phi ~ 0.64
# that collapses the genome-wide Empirical Bayes prior.
@inline function _is_polymorphic_event(ev::GLMEventData)::Bool
   sum_k = 0.0
   sum_r = 0.0
   sum_n = 0.0
   min_psi = Inf
   max_psi = -Inf
   @inbounds for j in 1:length(ev.psi)
      p = ev.psi[j]
      n = ev.reads[j]
      if !isfinite(p) || !isfinite(n) || n <= 0.0
         continue
      end
      kj = p * n
      sum_k += kj
      sum_r += n - kj
      sum_n += n
      p < min_psi && (min_psi = p)
      p > max_psi && (max_psi = p)
   end
   (sum_n <= 0.0 || sum_k < 1.0 || sum_r < 1.0) && return false
   pooled_psi = sum_k / sum_n
   (pooled_psi <= 0.01 || pooled_psi >= 0.99) && return false
   return (max_psi - min_psi) >= 1e-4
end

@inline function _is_utr_event(ev::GLMEventData)::Bool
   return ev.etype == "TS" || ev.etype == "TE"
end

# Run full 3-pass Beta-Binomial GLM + Empirical Bayes dispersion shrinkage across a cohort
function fit_glm_cohort(events::Vector{GLMEventData},
                        design::DesignSpec;
                        shrink_dispersion::Bool=true,
                        stratify_eb::Bool=true,
                        fixed_phi::Float64=-1.0,
                        ridge::Float64=1e-4)
   nev = length(events)
   results = Vector{GLMResult}(undef, nev)
   nev == 0 && return results, Float64[], (NaN, NaN, 0, nothing)

   # Check whether any event has residual degrees of freedom (m > p)
   has_replicates = false
   for ev in events
      mats = _build_event_matrices(design.X, ev.valid_idx, design.test_cols)
      if mats !== nothing && length(ev.valid_idx) > size(mats[1], 2)
         has_replicates = true
         break
      end
   end

   if fixed_phi > 0.0
      # User-specified fixed dispersion
      th_fixed = clamp(log((1.0 - fixed_phi) / fixed_phi), _THETA_MIN, _THETA_MAX)
      Threads.@threads for i in 1:nev
         results[i] = _test_single_event(events[i], design, th_fixed; ridge=ridge)
      end
      pvals = [res.pvalue for res in results]
      return results, bh_fdr(pvals), (th_fixed, 0.0, 0, nothing)
   end

   if !has_replicates
      println(stderr, "Warning: No residual degrees of freedom (unreplicated design); using default dispersion phi = 0.01.")
      th_def = log(99.0) # phi = 0.01
      Threads.@threads for i in 1:nev
         results[i] = _test_single_event(events[i], design, th_def; ridge=ridge)
      end
      pvals = [res.pvalue for res in results]
      return results, bh_fdr(pvals), (th_def, 0.0, 0, nothing)
   end

   # Pass 2: Raw Cox-Reid dispersion estimation strictly on polymorphic events
   # (skipping constitutive/boundary events whose flat likelihood creates artificial phi ~ 0.64 peaks).
   raw_thetas  = fill(NaN, nev)
   is_interior = fill(false, nev)
   samp_vars   = fill(NaN, nev)
   res_dfs     = zeros(Int, nev)
   is_poly     = fill(false, nev)

   Threads.@threads for i in 1:nev
      ev = events[i]
      _is_polymorphic_event(ev) || continue
      is_poly[i] = true
      mats = _build_event_matrices(design.X, ev.valid_idx, design.test_cols)
      mats === nothing && continue
      X_full = mats[1]
      m = length(ev.valid_idx)
      d_res = m - size(X_full, 2)
      d_res > 0 || continue

      k = Vector{Float64}(undef, m)
      r = Vector{Float64}(undef, m)
      n = ev.reads
      @inbounds for j in 1:m
         kj = ev.psi[j] * n[j]
         k[j] = kj
         r[j] = n[j] - kj
      end
      th, intr, sv = estimate_raw_theta(k, r, n, X_full; ridge=ridge)
      raw_thetas[i]  = th
      is_interior[i] = intr
      samp_vars[i]   = sv
      res_dfs[i]     = d_res
   end

   combined_prior = fit_eb_dispersion_prior(raw_thetas, is_interior, samp_vars, res_dfs)

   sp_prior = combined_prior
   utr_prior = nothing
   if stratify_eb
      sp_th  = Float64[]; sp_int = Bool[]; sp_sv = Float64[]; sp_df = Int[]; n_sp_int = 0
      u_th   = Float64[]; u_int  = Bool[]; u_sv  = Float64[]; u_df  = Int[]; n_u_int  = 0
      @inbounds for i in 1:nev
         is_poly[i] || continue
         if _is_utr_event(events[i])
            push!(u_th, raw_thetas[i]); push!(u_int, is_interior[i])
            push!(u_sv, samp_vars[i]);  push!(u_df, res_dfs[i])
            if is_interior[i] && isfinite(raw_thetas[i])
               n_u_int += 1
            end
         else
            push!(sp_th, raw_thetas[i]); push!(sp_int, is_interior[i])
            push!(sp_sv, samp_vars[i]);  push!(sp_df, res_dfs[i])
            if is_interior[i] && isfinite(raw_thetas[i])
               n_sp_int += 1
            end
         end
      end
      if n_sp_int >= 20
         sp_prior = fit_eb_dispersion_prior(sp_th, sp_int, sp_sv, sp_df)
      end
      if n_u_int >= 20
         utr_prior = fit_eb_dispersion_prior(u_th, u_int, u_sv, u_df)
      end
   end

   prior_mu, prior_var, n_prior = sp_prior
   # Effective prior degrees of freedom from Var(log chi2_d0) ≈ 2 / d0
   df_prior_sp  = shrink_dispersion ? clamp(2.0 / prior_var, 1.0, 30.0) : 0.0
   df_prior_utr = utr_prior !== nothing ? (shrink_dispersion ? clamp(2.0 / utr_prior[2], 1.0, 30.0) : 0.0) : df_prior_sp

   # Pass 3: EB MAP dispersion shrinkage + LRT/F-test per event in parallel
   Threads.@threads for i in 1:nev
      ev = events[i]
      mats = _build_event_matrices(design.X, ev.valid_idx, design.test_cols)
      if mats === nothing
         results[i] = _EMPTY_GLM_RESULT
         continue
      end
      X_full = mats[1]
      m = length(ev.valid_idx)
      d_res = m - size(X_full, 2)
      has_df = d_res > 0

      if _is_utr_event(ev) && utr_prior !== nothing
         pmu_s, pvar_s, df_p_s = utr_prior[1], utr_prior[2], df_prior_utr
      else
         pmu_s, pvar_s, df_p_s = prior_mu, prior_var, df_prior_sp
      end

      th_use = pmu_s
      if is_poly[i] && has_df
         if shrink_dispersion
            k = Vector{Float64}(undef, m)
            r = Vector{Float64}(undef, m)
            n = ev.reads
            @inbounds for j in 1:m
               kj = ev.psi[j] * n[j]
               k[j] = kj
               r[j] = n[j] - kj
            end
            th_use = estimate_eb_theta(k, r, n, X_full, pmu_s, pvar_s; ridge=ridge)
         elseif isfinite(raw_thetas[i])
            th_use = raw_thetas[i]
         end
      end

      df_denom = (d_res > 0 || df_p_s > 0.0) ? Float64(max(d_res, 0)) + df_p_s : Inf
      results[i] = _test_single_event(ev, design, th_use; df_denom=df_denom, ridge=ridge)
   end

   pvals = [res.pvalue for res in results]
   fdr = bh_fdr(pvals)
   return results, fdr, (prior_mu, prior_var, n_prior, utr_prior)
end

@inline function _fmt_num(x::Float64, digits::Int)::String
   isfinite(x) ? string(round(x, digits=digits)) : "NA"
end

@inline function _fmt_pval(x::Float64)::String
   if !isfinite(x)
      return "NA"
   elseif x == 0.0
      return "0.0"
   elseif x < 1e-4
      return @sprintf("%.4e", x)
   else
      return @sprintf("%.6f", x)
   end
end

function output_diff_glm(outfile::String,
                         events::Vector{GLMEventData},
                         results::Vector{GLMResult},
                         fdr::Vector{Float64},
                         design::DesignSpec;
                         sig::Int=5)
   io = open(outfile, "w")
   stream = ZlibDeflateOutputStream(io)
   if design.is_two_group
      write(stream, "Gene\tNode\tCoord\tStrand\tType\tComplexity\tEntropy\t" *
                    "SampleCnt_A\tSampleCnt_B\tMeanReads_A\tMeanReads_B\t" *
                    "Psi_A\tPsi_B\tDeltaPsi\tLogOddsRatio\tSE\tDispersion\tLRT\tPvalue\tFDR\n")
      @inbounds for i in 1:length(events)
         res = results[i]
         res.tested || continue
         ev = events[i]
         write(stream, ev.gene, '\t', ev.node, '\t', ev.coord, '\t', ev.strand, '\t', ev.etype, '\t',
                       COMPLEX_CHAR, string(ev.complex), '\t', string(round(ev.entropy, digits=3)), '\t',
                       string(ev.samp_a), '\t', string(ev.samp_b), '\t',
                       string(round(ev.mean_reads_a, digits=1)), '\t',
                       string(round(ev.mean_reads_b, digits=1)), '\t',
                       _fmt_num(res.psi_a, sig), '\t',
                       _fmt_num(res.psi_b, sig), '\t',
                       _fmt_num(res.delta_psi, sig), '\t',
                       _fmt_num(res.log_odds, sig), '\t',
                       _fmt_num(res.se, sig), '\t',
                       _fmt_pval(res.dispersion), '\t',
                       _fmt_num(res.lrt, 4), '\t',
                       _fmt_pval(res.pvalue), '\t',
                       _fmt_pval(fdr[i]), '\n')
      end
   else
      write(stream, "Gene\tNode\tCoord\tStrand\tType\tComplexity\tEntropy\t" *
                    "SampleCnt\tMeanReads\tMeanPsi\tDeltaPsi\tLogOddsRatio\tSE\t" *
                    "Dispersion\tLRT\tDF\tPvalue\tFDR\n")
      @inbounds for i in 1:length(events)
         res = results[i]
         res.tested || continue
         ev = events[i]
         write(stream, ev.gene, '\t', ev.node, '\t', ev.coord, '\t', ev.strand, '\t', ev.etype, '\t',
                       COMPLEX_CHAR, string(ev.complex), '\t', string(round(ev.entropy, digits=3)), '\t',
                       string(ev.samp_a), '\t',
                       string(round(ev.mean_reads_a, digits=1)), '\t',
                       _fmt_num(res.psi_a, sig), '\t',
                       _fmt_num(res.delta_psi, sig), '\t',
                       _fmt_num(res.log_odds, sig), '\t',
                       _fmt_num(res.se, sig), '\t',
                       _fmt_pval(res.dispersion), '\t',
                       _fmt_num(res.lrt, 4), '\t',
                       string(res.df), '\t',
                       _fmt_pval(res.pvalue), '\t',
                       _fmt_pval(fdr[i]), '\n')
      end
   end
   close(stream)
   close(io)
end

# High-level entry point called by bin/whippet-delta.jl
function process_psi_files_glm(outfile::String,
                               design::DesignSpec;
                               min_samp::Int=1,
                               min_reads::Int=5,
                               shrink_dispersion::Bool=true,
                               stratify_eb::Bool=true,
                               fixed_phi::Float64=-1.0,
                               ridge::Float64=1e-4)
   events = read_psi_cohort(design; min_reads=min_reads, min_samp=min_samp)
   println(stderr, "Loaded $(length(events)) testable events across $(length(design.files)) samples (using $(Threads.nthreads()) thread(s))...")
   results, fdr, (pmu, pvar, np, utr_prior) = fit_glm_cohort(events, design;
                                                             shrink_dispersion=shrink_dispersion,
                                                             stratify_eb=stratify_eb,
                                                             fixed_phi=fixed_phi,
                                                             ridge=ridge)
   if fixed_phi <= 0.0 && isfinite(pmu)
      prior_phi = 1.0 / (1.0 + exp(pmu))
      d0_sp = shrink_dispersion ? clamp(2.0 / pvar, 1.0, 30.0) : 0.0
      if utr_prior !== nothing
         umu, uvar, unp = utr_prior
         uphi = 1.0 / (1.0 + exp(umu))
         d0_u = shrink_dispersion ? clamp(2.0 / uvar, 1.0, 30.0) : 0.0
         println(stderr, "Empirical Bayes dispersion prior [Spliced] (from $np events): median phi = $(round(prior_phi, sigdigits=3)), sd(log M) = $(round(sqrt(pvar), digits=3)), d0 = $(round(d0_sp, digits=2))")
         println(stderr, "Empirical Bayes dispersion prior [Tandem UTR] (from $unp events): median phi = $(round(uphi, sigdigits=3)), sd(log M) = $(round(sqrt(uvar), digits=3)), d0 = $(round(d0_u, digits=2))")
      else
         println(stderr, "Empirical Bayes dispersion prior (from $np events): median phi = $(round(prior_phi, sigdigits=3)), sd(log M) = $(round(sqrt(pvar), digits=3)), d0 = $(round(d0_sp, digits=2))")
      end
   end
   n_tested = count(r -> r.tested, results)
   n_sig = count(i -> results[i].tested && isfinite(fdr[i]) && fdr[i] <= 0.05, 1:length(results))
   println(stderr, "Tested $n_tested events ($n_sig with FDR <= 0.05).")
   output_diff_glm(outfile, events, results, fdr, design)
   return events, results, fdr
end
