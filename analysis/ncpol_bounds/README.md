# Bounded ncpol2sdpa diagnostics

All entries in the table below are numerical diagnostics, **certified_exact=false**. No threshold or exact-level assignment follows from a small gap. Two separated proposals have since been converted to [exact rational witnesses](rational_README.md) and [certified failure intervals](INTERVAL_PROOF.md). Those proofs have separate solver-free replays.

|family|alpha|words|status|primal max|dual max|quantum reference|primal minus Q|
|---|---:|---:|---|---:|---:|---:|---:|
|symmetric|0|9|optimal_inaccurate|2.82842712211|2.8284271273|2.82842712475|-2.639e-09|
|symmetric|0|13|optimal_inaccurate|2.8284271232|2.82842713068|2.82842712475|-1.551e-09|
|symmetric|1/4|9|optimal_inaccurate|2.90449094405|2.904490961|2.90449095075|-6.697e-09|
|symmetric|1/4|13|optimal_inaccurate|2.90449091469|2.90449097147|2.90449095075|-3.606e-08|
|symmetric|1/2|9|optimal_inaccurate|3.13693729401|3.13693732745|3.13693731131|-1.73e-08|
|symmetric|1/2|13|optimal_inaccurate|3.13693721199|3.13693735986|3.13693731131|-9.932e-08|
|symmetric|3/4|9|optimal|3.52013366873|3.52013366874|3.51907732176|0.001056|
|symmetric|3/4|13|optimal_inaccurate|3.51907730452|3.5190773309|3.51907732176|-1.725e-08|
|symmetric|19/20|9|optimal_inaccurate|3.90051761245|3.90051761245|3.90016389937|0.0003537|
|symmetric|19/20|13|optimal_inaccurate|3.90032966397|3.90032968159|3.90016389937|0.0001658|
|one_sided|1|5|optimal|3.82842712461|3.82842712496|3.16227766017|0.6661|
|one_sided|1|9|optimal_inaccurate|3.16227765704|3.16227766493|3.16227766017|-3.13e-09|

Full word lists, versions, solver settings and residuals are in bounds.json.
Only universal involutions and cross-party commutation enter the moments.
The two objective values are reconstructed separately from the CVXPY primal and PSD dual; ncpol2sdpa's solve wrapper would duplicate one value.
Independent 1+AB assembly crosscheck at alpha=1/2: absolute difference 1.101e-10.
Total batch wall time: 0.425 s.
