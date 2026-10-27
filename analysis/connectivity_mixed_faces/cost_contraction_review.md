# Review of contraction from a local quantum optimum

Anubhav Chaturvedi

28 October 2026. The contraction lemma is correct. It proves that every path component of the full 1+AB support-gap domain has a representative with strict quantum advantage over the local bound. It does not prove that the strict-advantage subset is connected, or that its components remain distinct after paths through local-optimum directions are allowed.

Fix a CHSH sector and its eight losing-event probabilities p_i. Write

    m_T(c)=min_(q in T) sum_i c_i p_i(q),    T=Q or Q1AB.

Every excluded direction has c_i>0. Positive scaling permits min_i c_i=1. Suppose its quantum optimum equals the local optimum, m_Q(c)=1, and choose a 1+AB witness x with c.p(x)<1. Set

    d=c-1 >=0,    c(t)=1+t d,    delta=1-c.p(x)>0.

For all t in [0,1], all costs are positive and their minimum is one. A deterministic vertex corresponding to any index with d_i=0 has loss one, while every other deterministic loss is at least one. Thus the local optimum stays exactly one. The minimum-index set is unchanged for **t>0**; at t=0 all indices are minima. This harmless exception should be stated correctly.

Probability positivity makes m_Q(c(t)) nondecreasing. It is continuous because Q is compact and the objective depends continuously and linearly on t. The CHSH value gives

    m_Q(c(0))=2-sqrt(2)<1,    m_Q(c(1))=1,
    m_Q(c(t))<=1.

Consequently there is a first t_Q in (0,1], and

    {t in [0,1]:m_Q(c(t))=1}=[t_Q,1].

This includes the possible endpoint case t_Q=1. No strict monotonicity or differentiability is needed.

The original witness works uniformly on the whole contraction:

    m_Q1AB(c(t)) <= c(t).p(x) <= c(1).p(x)=1-delta,

because d_i>=0 and p_i(x)>=0. Hence for every t in [t_Q,1], the support gap in loss convention is at least delta. Continuity at t_Q supplies t_- with 0<t_-<t_Q such that

    m_Q(c(t))>1-delta/2  for t in [t_-,t_Q].

On that interval the same witness gives a gap greater than delta/2. Also m_Q(c(t))<1 for t<t_Q, by definition of the first threshold. Thus c(t_-) has strict quantum advantage and the entire segment from c(1) to c(t_-) remains excluded.

To recover the chosen normalized Bell chart, replace c(t) by c(t)/sum_i c_i(t). This continuous positive scaling preserves both strict inequalities. All costs remain positive, so the path stays in the same unique-PR sector. The H16 quotient maps it to a path with the same endpoint properties. Already strict-advantage excluded directions require no contraction.

The previous deterministic perturbation lemma establishes open excluded sets with a unique local quantum maximizer. This contraction is compatible with that fact: those sets are not separate components consisting entirely of local-optimum directions. Paths are allowed to leave them while retaining the strict 1+AB gap. The remaining topology can therefore be studied using strict-advantage representatives, while retaining paths and adjacencies through the local-tie and local-only parts of the full gap domain.

No solver or new literature input is needed for this review. The only special quantum value used is the exact CHSH bound at the all-equal cost vector.
