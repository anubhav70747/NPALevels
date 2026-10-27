# Every failure component reaches strict quantum advantage

Anubhav Chaturvedi

Let K=Q and L=Q1+AB. Every component of the strict support-gap domain, before or after physical relabeling, contains a direction with strict quantum advantage over the local polytope. Thus components containing only local quantum optimizers do not form a separate unresolved class. The proof provides a path with one fixed supraquantum witness and a positive margin throughout.

## Cost coordinates and a quantitative path

Every excluded objective has a unique PR maximizing sector. In its event-cost chart write

    b.R-b.q = sum_i c_i p_i(q),    c_i>0.

Rescale to min_i c_i=1, and set a_i=c_i-1>=0. A deterministic neighbor whose cost is one attains the local minimum of the cost functional. The other deterministic vertices lose three events, so cannot improve that minimum. Write

    m_T(c)=min_{q in T} sum_i c_i p_i(q).

Support inexactness means m_K(c)>m_L(c); quantum advantage over the local value means m_K(c)<1. If advantage is already strict, no deformation is needed.

Otherwise m_K(c)=1. Choose a strict witness x in L and put

    delta = 1-sum_i c_i p_i(x)>0,
    c(t) = 1+t a,       0<=t<=1.                       (1)

All costs remain positive. For every t>0 their minimum is one and their minimizing index set is unchanged. Since event probabilities are nonnegative,

    sum_i c_i(t)p_i(x) <= 1-delta.                       (2)

The function m_K(c(t)) is continuous and nondecreasing, and is at most one because the same deterministic vertex attains one. At t=0 the costs are uniform. The identity sum_i p_i=2-CHSH/2 and the exact quantum CHSH maximum give

    m_K(c(0))=2-sqrt(2)=m_L(c(0)).                       (3)

Moreover (2) at t=0 implies m_K(c(0))<=1-delta. At t=1 the quantum minimum is one. Hence continuity supplies t_delta in (0,1) with

    m_K(c(t_delta))=1-delta/2.                           (4)

Take, for example, the first such t. Monotonicity and (2) now give, for EVERY t in [t_delta,1],

    m_K(c(t))-sum_i c_i(t)p_i(x) >= delta/2>0.            (5)

The reversed segment from t=1 to t=t_delta is the required path. Its endpoint has quantum advantage delta/2 and still has a certified witness margin at least delta/2. This is stronger than a continuity assertion at the first quantum/local transition. It does not claim that the complete segment to uniform CHSH stays excluded; that endpoint is admitted.

To return to the normalized cost simplex divide by sum_i c_i(t). This is a continuous positive rescaling, and sum_i c_i(t)<=sum_i c_i. Thus the normalized margin is at least delta/(2 sum_i c_i). The path never changes its unique PR sector or crosses a zero-cost face. The construction is covariant under the sixteen residual event permutations and therefore descends to the physical symmetry quotient.

## Scope and subsequent routing

The choice of t_delta uses the exact quantum minimum along a one-dimensional cost path. It is an existence construction; no numerical root or universal efficient algorithm is asserted. For algebraic input and witness it can in principle be selected by the existing exact semialgebraic quantum predicate.

The argument also holds for each fixed standard NPA level containing 1+AB: event probabilities are nonnegative and that level has the same exact uniform-cost CHSH bound. The selected path and endpoint may depend on the chosen level and witness. No simultaneous all-level path or monotonicity of the minimum exact depth is claimed.

At the resulting strict-advantage objective, a pure real-qubit strategy attains the quantum maximum by the compact atom model. Its behavior is nonlocal, so neither local measurement pair commutes. It therefore supplies a finite noncommuting stationary contact for [the canonical-source route](normal_routes.md). This step needs no uniqueness assumption and no assertion that all nonlocal exposed contacts have positive probabilities.

Combining the two routes reduces the universal Motzkin connection problem to the excluded parts of the two canonical C1/C2 source families, with their boundaries retained. The contraction alone does not connect those remaining families to the Motzkin collar.
