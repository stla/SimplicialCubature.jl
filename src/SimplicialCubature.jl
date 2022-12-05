module SimplicialCubature

export CanonicalSimplex
export integrateOnSimplex
export integratePolynomialOnSimplex

import LinearAlgebra
import TypedPolynomials

function SimplexVolume(S)
  n, _ = size(S)
  v = S[:,n+1]
  V = Array{Float64}(undef, n, 0)
  for i = 1:n
    V = hcat(V, S[:,i] - v)
  end
  return abs(LinearAlgebra.det(V)) / factorial(n)
end

function vectorOfVectorsToMatrix(V)
  return hcat(V...)
end

"""
    CanonicalSimplex(n)

Canonical n-dimensional simplex.

# Argument
- `n`: positive integer
"""
function CanonicalSimplex(n)
  S = hcat(fill(0.0, n, 1), LinearAlgebra.diagm(fill(1.0, n)))
  return map(i -> S[:,i], 1:(n+1))
end

"""
    integrateOnSimplex(f, S; dim, maxEvals, absError, tol, rule, info, fkwargs...)

Integration of a function over one or more simplices.

# Arguments
- `f`: function to be integrated; must return a real scalar value or a real vector
- `S`: simplex or vector of simplices; a simplex is given by n+1 vectors of dimension n
- `dim`: number of components of `f`
- `maxEvals`: maximum number of calls to `f`
- `absError`: requested absolute error
- `tol`: requested relative error
- `rule`: integration rule, an integer between 1 and 4; a 2*rule+1 degree integration rule is used
- `info`: Boolean, whether to print more info
- `fkwargs`: keywords arguments of `f`
"""
function integrateOnSimplex(
  f::Function,
  S::Union{Vector{Vector{T}},Vector{Vector{Vector{T}}}};
  dim = 1,
  maxEvals = 10000,
  absError = 0.0,
  tol = 1.0e-5,
  rule = 3,
  info = false,
  fkwargs...,
) where T <: Real
  t = eltype(eltype(S))
  if t == eltype(t)
    S = [S]
  end
  nS = length(S)
  m = length(S[1])
  n = length(S[1][1])
  if m != n + 1
    error("Invalid simplex")
  end
  Simplices = Vector{Array{Float64, 2}}(undef, nS)
  for i in 1:nS
    Si = S[i]
    if length(Si) != m
      error("Invalid simplex found in `S`.")
    end
    for j in 1:m
      if length(Si[j]) != n
        error("Invalid simplex found in `S`.")
      end
    end
    Simplices[i] = vectorOfVectorsToMatrix(convert(Vector{Vector{Float64}}, Si))
  end
  function fNew(x)
    return [f(x; fkwargs...)]
  end
  a = adsimp(n, Simplices, dim, fNew, maxEvals, absError, tol, rule, info)
  local result
  if info
    result = (
      integral = a.VL,
      estAbsError = a.AE,
      functionEvaluations = a.NV,
      returnCode = a.FL,
      subsimplices = a.VRTS,
      subsimplicesIntegral = a.VLS,
      subsimplicesAbsError = a.AES,
      subsimplicesVolume = a.VOL,
      message = adsimp_message(a.FL),
    )
  else
    result = (
      integral = a.VL,
      estAbsError = a.AE,
      functionEvaluations = a.NV,
      returnCode = a.FL,
      message = adsimp_message(a.FL),
    )
  end
  return result
end

function adsimp_message(rcode)
  local msg
  if rcode == 0 
    msg = "OK"
  elseif rcode == 1 
    msg = "error: maxEvals exceeded - too many function evaluations"
  elseif rcode == 2
    msg = "error: integRule < 0 or integRule > 4"
  elseif rcode == 3
    msg = "error: dimension n of the space < 2"
  elseif rcode == 4
    msg = "error: fDim = dimension of f < 1"
  elseif rcode == 5
    msg = "error: absError < 0 and tol < 0"
  else
    msg = "error: unknown return code = $rcode"
  end 
  return msg
end

function adsimp(ND, VRTS, NF, F, MXFS, EA, ER, KEY, partitionInfo)
  if KEY == 0
    KEY = 3
  end
  SBS = length(VRTS)
  b = SMPCHC(ND, NF, MXFS, EA, ER, SBS, KEY)
  if b.FL != 0
    return (
      VL = fill(NaN, NF),
      AE = fill(Inf, NF),
      NV = 0,
      FL = b.VL,
      VLS = NaN,
      AES = NaN,
      VOL = NaN,
    )
  end
  return SMPSAD(
    ND,
    NF,
    F,
    MXFS,
    EA,
    ER,
    KEY,
    b.RULCLS,
    SBS,
    VRTS,
    partitionInfo,
  )
end

function SMPCHC(ND, NF, MXFS, EA, ER, SBS, KEY)
  FL = 0
  if (KEY < 0) || (KEY > 4)
    FL = 2
  end
  if ND < 2
    FL = 3
  end
  if NF < 1
    FL = 4
  end
  if (EA < 0) && (ER < 0)
    FL = 5
  end
  local RULCLS
  if FL == 0
    if KEY == 0
      RULCLS = (ND + 4) * (ND + 3) * (ND + 2) / 6 + (ND + 2) * (ND + 1)
    elseif KEY == 1
      RULCLS = 2 * ND + 3
    elseif KEY == 2
      RULCLS = (ND + 3) * (ND + 2) / 2 + 2 * (ND + 1)
    elseif KEY == 3
      RULCLS = (ND + 4) * (ND + 3) * (ND + 2) / 6 + (ND + 2) * (ND + 1)
    elseif KEY == 4
      RULCLS =
        (ND + 5) * (ND + 4) * (ND + 3) * (ND + 2) / 24 +
        5 * (ND + 2) * (ND + 1) / 2
    end
    if MXFS < SBS * RULCLS
      FL = 1
    end
  else
    RULCLS = 0
  end
  return (FL = FL, RULCLS = RULCLS)
end

function SMPSAD(ND, NF, F, MXFS, EA, ER, KEY, RCLS, SBS, VRTS, partitionInfo)
  NV = 0
  DFCOST = 1 + 2 * ND * (ND + 1)
  VL = fill(0.0, NF)
  AE = VL
  a = SMPRMS(ND, KEY)
  FCT = factorial(ND)
  VLS = fill(0.0, NF, SBS) # VLS[i,j] = estimated integral of F[i] on simplex VRTS[,,j]
  AES = fill(0.0, NF, SBS) # AES[i,j] = estimated abs. err. in integral of F[i] on simplex VRTS[,,j]
  VOL = fill(0.0, SBS)
  for K = 1:SBS
    VOL[K] =
      abs(
        LinearAlgebra.det(VRTS[K][:, 1:ND] - repeat(VRTS[K][:, ND+1], 1, ND)),
      ) / FCT
    b = SMPRUL(ND, VRTS[K], VOL[K], NF, F, a.G, a.W, a.PTS)
    AES[:, K] = b.RGNERR
    VLS[:, K] = b.BASVAL
    VL = VL + VLS[:, K]
    AE = AE + AES[:, K]
    NV = NV + RCLS
  end
  FL = convert(Int64, maximum(AE .> maximum([EA, (ER .* abs.(VL))...])))
  while ((FL > 0) && (NV + DFCOST + 4 * RCLS <= MXFS))
    ID = argmax(map(maximum, eachcol(AES)))
    VL = VL - VLS[:, ID]
    AE = AE - AES[:, ID]
    (VRTS, NEW) = SMPDFS(ND, NF, F, ID, SBS, VRTS)
    VI = VOL[ID] / NEW
    VOL = [VOL..., fill(0.0, NEW - 1)...]
    VLS = hcat(VLS, fill(0.0, NF, NEW - 1))
    AES = hcat(AES, fill(0.0, NF, NEW - 1))
    for K in [ID, collect((SBS+1):(SBS+NEW-1))...]
      VOL[K] = VI
      d = SMPRUL(ND, VRTS[K], VI, NF, F, a.G, a.W, a.PTS)
      VLS[:, K] = d.BASVAL
      AES[:, K] = d.RGNERR
      VL = VL + VLS[:, K]
      AE = AE + AES[:, K]
      NV = NV + RCLS
    end
    NV = NV + DFCOST
    SBS = SBS + NEW - 1
    FL = convert(Int64, maximum(AE .> maximum([EA, (ER .* abs.(VL))...])))
  end
  if SBS > 1
    VL = map(sum, eachrow(VLS))
    AE = map(sum, eachrow(AES))
  end
  if NF == 1
    VL = VL[1]
    AE = AE[1]
  end
  local result
  if partitionInfo
    result = (
      VL = VL,
      AE = AE,
      NV = NV,
      FL = FL,
      VRTS = VRTS,
      VLS = VLS,
      AES = AES,
      VOL = VOL,
    )
  else
    result = (VL = VL, AE = AE, NV = NV, FL = FL)
  end
  return result
end

function SMPDFS(ND, NF, F, TOP, SBS, VRTS)
  CUTTF = 2
  CUTTB = 8
  IS = 1
  JS = 2
  DFMX = 0
  EMX = 0
  V = VRTS[TOP]
  _, n = size(V)
  CN = map(sum, eachrow(V)) ./ n
  FC = F(CN)
  DFMD = sum(abs.(FC))
  FRTHDF = fill(0.0, ND, ND + 1)
  local IE, JE, IT, JT, DFNX
  for I = 1:ND
    for J = (I+1):(ND+1)
      H = 2 * (V[:, I] - V[:, J]) ./ (5 * (ND + 1))
      EWD = sum(abs.(H))
      if EWD >= EMX
        IE = I
        JE = J
        EMX = EWD
      end
      DFR = sum(
        abs.(
          F(CN - 2 * H) + F(CN + 2 * H) + 6 * FC - 4 * (F(CN - H) + F(CN + H))
        ),
      )
      if (DFMD + DFR / 8) == DFMD
        DFR = 0
      end
      DFR = DFR * EWD
      if DFR >= DFMX
        IT = IS
        JT = JS
        DFNX = DFMX
        IS = I
        JS = J
        DFMX = DFR
      else
        if DFR >= DFNX # ??? DFNX not defined !
          IT = I
          JT = J
          DFNX = DFR
        end
      end
      FRTHDF[I, J] = DFR
    end
  end
  local NEW, LS
  if DFNX > DFMX / CUTTF
    NEW = 4
  else
    NEW = 3
    if (DFMX == 0)
      IS = IE
      JS = JE
    else
      DFSMX = 0
      for L = 1:(ND+1)
        if (L != IS) && (L != JS)
          IT = minimum([L, IS, JS])
          JT = maximum([L, IS, JS])
          LT = IS + JS + L - IT - JT
          DFR = FRTHDF[IT, LT] + FRTHDF[LT, JT]
          if DFR >= DFSMX
            DFSMX = DFR
            LS = L
          end
        end
      end
      DIFIL = FRTHDF[min(IS, LS), max(IS, LS)]
      DIFLJ = FRTHDF[min(JS, LS), max(JS, LS)]
      DFNX = DIFIL + DIFLJ - min(DIFIL, DIFLJ)
      if (DFMX / CUTTB < DFNX) && (DIFIL > DIFLJ)
        IT = IS
        IS = JS
        JS = IT
      end
    end
  end
  VV = fill(V, NEW - 1)
  VRTS = vcat(VRTS, VV)
  VTI = V[:, IS]
  VTJ = V[:, JS]
  if NEW == 4
    VRTS[TOP][:, JS] = (VTI + VTJ) / 2
    VRTS[SBS+1][:, IS] = VTI
    VRTS[SBS+1][:, JS] = (VTI + VTJ) / 2
    VRTS[SBS+2][:, IS] = (VTI + VTJ) / 2
    VRTS[SBS+2][:, JS] = VTJ
    VRTS[SBS+3][:, IS] = (VTI + VTJ) / 2
    VRTS[SBS+3][:, JS] = VTJ
    VTI = VRTS[TOP][:, IT]
    VTJ = VRTS[TOP][:, JT]
    VRTS[TOP][:, JT] = (VTI + VTJ) / 2
    VRTS[SBS+1][:, IT] = (VTI + VTJ) / 2
    VRTS[SBS+1][:, JT] = VTJ
    VTI = VRTS[SBS+2][:, IT]
    VTJ = VRTS[SBS+2][:, JT]
    VRTS[SBS+2][:, JT] = (VTI + VTJ) / 2
    VRTS[SBS+3][:, IT] = (VTI + VTJ) / 2
    VRTS[SBS+3][:, JT] = VTJ
  else
    VRTS[TOP][:, JS] = (2 * VTI + VTJ) / 3
    VRTS[SBS+1][:, IS] = (2 * VTI + VTJ) / 3
    if DFMX / CUTTF < DFNX
      VRTS[SBS+1][:, JS] = VTJ
      VRTS[SBS+2][:, IS] = (2 * VTI + VTJ) / 3
      VRTS[SBS+2][:, JS] = VTJ
      VTJ = VRTS[SBS+1][:, JS]
      VTL = VRTS[SBS+1][:, LS]
      VRTS[SBS+1][:, LS] = (VTJ + VTL) / 2
      VRTS[SBS+2][:, JS] = (VTJ + VTL) / 2
      VRTS[SBS+2][:, LS] = VTL
    else
      VRTS[SBS+1][:, JS] = (VTI + 2 * VTJ) / 3
      VRTS[SBS+2][:, IS] = (VTI + 2 * VTJ) / 3
      VRTS[SBS+2][:, JS] = VTJ
    end
  end
  return (VRTS = VRTS, NEW = NEW)
end

function SMPRUL(ND, VRTS, VOL, NF, F, G, W, PTS)
  RTMN = 0.1
  SMALL = 1.0e-12
  ERRCOF = 8
  WTS, RLS = size(W)
  RULE = fill(0.0, NF, RLS)
  for K = 1:WTS
    if PTS[K] > 0
      RULE =
        RULE .+
        VOL *
        SMPSMS(ND, VRTS, NF, F, G[:, K]) *
        LinearAlgebra.transpose(W[K, :])
    end
  end
  BASVAL = fill(0.0, NF)
  RGNERR = fill(0.0, NF)
  local NMCP
  for I = 1:NF
    BASVAL[I] = RULE[I, 1]
    NMBS = abs(BASVAL[I])
    RT = RTMN
    K = RLS
    while K >= 3
      NMRL = max(abs(RULE[I, K]), abs(RULE[I, K-1]))
      if (NMRL > SMALL * NMBS) && (K < RLS)
        RT = max(NMRL / NMCP, RT)
      end
      RGNERR[I] = max(NMRL, RGNERR[I])
      NMCP = NMRL
      K = K - 2
    end
    if (RT < 1) && (RLS > 3)
      RGNERR[I] = RT * NMCP
    end
    RGNERR[I] = max(ERRCOF * RGNERR[I], SMALL * NMBS)
  end
  return (BASVAL = BASVAL, RGNERR = RGNERR)
end

function SMPSMS(N, VERTEX, NF, F, G)
  SYMSMS = fill(0.0, NF)
  G = sort(G, rev = true)
  pr = true
  local LX
  while pr
    SYMSMS = SYMSMS + F(VERTEX * G)
    pr = false
    for I = 2:(N+1)
      GI = G[I]
      if G[I-1] > GI
        IX = I - 1
        for L = 1:div(IX + 1, 2) # !!!!!!!!!!
          GL = G[L]
          if GL <= GI
            IX = IX - 1
          end
          G[L] = G[I-L]
          G[I-L] = GL
          if G[L] > GI
            LX = L
          end
        end
        if G[IX] <= GI
          IX = LX
        end
        G[I] = G[IX]
        G[IX] = GI
        pr = true
        break
      end
    end
  end
  return repeat(SYMSMS, 1, 1)
end

function SMPRMS(N, KEY)
  local RLS, GMS, WTS
  if KEY == 1
    RLS = 3
    GMS = 2
    WTS = 3
  elseif KEY == 2
    RLS = 5
    GMS = 4
    WTS = 6
  elseif KEY == 3
    RLS = 7
    GMS = 7
    WTS = 11
  elseif KEY == 4
    RLS = 7
    GMS = 12
    WTS = 21
    if N == 2
      GMS = 11
      WTS = 20
    end
  end
  W = fill(0.0, WTS, RLS)
  PTS = fill(0.0, WTS)
  G = fill(0.0, N + 1, WTS)
  NP = N + 1
  N2 = NP * (N + 2)
  N4 = N2 * (N + 3) * (N + 4)
  N6 = N4 * (N + 5) * (N + 6)
  N8 = N6 * (N + 7) * (N + 8)
  G[:, 1] .= 1 / NP
  PTS[1] = 1
  R1 = (N + 4 - sqrt(15)) / (N * N + 8 * N + 1)
  S1 = 1 - N * R1
  L1 = S1 - R1
  G[1, GMS+1] = S1
  G[2:NP, GMS+1] .= R1
  PTS[GMS+1] = NP
  IW = RLS
  if KEY < 4
    W[1, IW] = 1
    IW = IW - 1
    W[GMS+1, IW] = 1 / NP
    IW = IW - 1
  end
  G[1, 2] = 3 / (N + 3)
  G[2:NP, 2] .= 1 / (N + 3)
  PTS[2] = NP
  W[2, IW] = (N + 3)^3 / (4 * N2 * (N + 3))
  if KEY > 1
    IW = IW - 1
    if N == 2
      L2 = 0.62054648267200632589046034361711
      L1 = -sqrt(1 / 2 - L2^2)
      R1 = (1 - L1) / 3
      S1 = 1 - 2 * R1
      G[1, GMS+1] = S1
      G[2:NP, GMS+1] .= R1
      PTS[GMS+1] = 3
      W[GMS+1, IW] = 1 / 6
      R2 = (1 - L2) / 3
      S2 = 1 - 2 * R2
      G[1, GMS+2] = S2
      G[2:NP, GMS+2] .= R2
      PTS[GMS+2] = 3
      W[GMS+2, IW] = 1 / 6
    else
      R2 = (N + 4 + sqrt(15)) / (N * N + 8 * N + 1)
      S2 = 1 - N * R2
      L2 = S2 - R2
      G[1, GMS+2] = S2
      G[2:NP, GMS+2] .= R2
      PTS[GMS+2] = NP
      W[GMS+2, IW] = (2 / (N + 3) - L1) / (N2 * (L2 - L1) * L2^2)
      W[GMS+1, IW] = (2 / (N + 3) - L2) / (N2 * (L1 - L2) * L1^2)
    end
    IW = IW - 1
    G[1, 3] = 5 / (N + 5)
    G[2:NP, 3] .= 1 / (N + 5)
    PTS[3] = NP
    G[1, 4] = 3 / (N + 5)
    G[2, 4] = 3 / (N + 5)
    G[3:NP, 4] .= 1 / (N + 5)
    PTS[4] = NP * N / 2
    W[2, IW] = -(N + 3)^5 / (16 * N4)
    W[3:4, IW] .= (N + 5)^5 / (16 * N4 * (N + 5))
  end
  if KEY > 2
    IW = IW - 1
    U1 = (N + 7 + 2 * sqrt(15)) / (N * N + 14 * N - 11)
    V1 = (1 - (N - 1) * U1) / 2
    D1 = V1 - U1
    G[1, GMS+3] = V1
    G[2, GMS+3] = V1
    G[3:NP, GMS+3] .= U1
    PTS[GMS+3] = ((N + 1) * N) / 2
    U2 = (N + 7 - 2 * sqrt(15)) / (N * N + 14 * N - 11)
    V2 = (1 - (N - 1) * U2) / 2
    D2 = V2 - U2
    G[1, GMS+4] = V2
    G[2, GMS+4] = V2
    G[3:NP, GMS+4] .= U2
    PTS[GMS+4] = ((N + 1) * N) / 2
    if N == 2
      W[GMS+3, IW] = (155 - sqrt(15)) / 1200
      W[GMS+4, IW] = (155 + sqrt(15)) / 1200
      W[1, IW] = 1 - 3 * (W[GMS+3, IW] + W[GMS+4, IW])
    elseif N == 3
      W[GMS+1, IW] = (2665 + 14 * sqrt(15)) / 37800
      W[GMS+2, IW] = (2665 - 14 * sqrt(15)) / 37800
      W[GMS+3, IW] = 2 * 15 / 567
      PTS[GMS+4] = 0
    else
      W[GMS+1, IW] =
        (2 * (27 - N) / (N + 5) - L2 * (13 - N)) / (L1^4 * (L1 - L2) * N4)
      W[GMS+2, IW] =
        (2 * (27 - N) / (N + 5) - L1 * (13 - N)) / (L2^4 * (L2 - L1) * N4)
      W[GMS+3, IW] = (2 / (N + 5) - D2) / (N4 * (D1 - D2) * D1^4)
      W[GMS+4, IW] = (2 / (N + 5) - D1) / (N4 * (D2 - D1) * D2^4)
    end
    IW = IW - 1
    G[1, 5] = 7 / (N + 7)
    G[2:NP, 5] .= 1 / (N + 7)
    PTS[5] = NP
    G[1, 6] = 5 / (N + 7)
    G[2, 6] = 3 / (N + 7)
    G[3:NP, 6] .= 1 / (N + 7)
    PTS[6] = NP * N
    G[1:3, 7] .= 3 / (N + 7)
    if NP > 3
      G[4:NP, 7] .= 1 / (N + 7)
    end
    PTS[7] = NP * N * (N - 1) / 6
    W[2, IW] = (N + 3)^7 / (2 * 64 * N4 * (N + 5))
    W[3:4, IW] .= -(N + 5)^7 / (64 * N6)
    W[5:7, IW] .= (N + 7)^7 / (64 * N6 * (N + 7))
  end
  if KEY == 4
    IW = IW - 1
    SG = 1 / (23328 * N6)
    U5 = -6^3 * SG * (52212 - N * (6353 + N * (1934 - N * 27)))
    U6 = 6^4 * SG * (7884 - N * (1541 - N * 9))
    U7 = -6^5 * SG * (8292 - N * (1139 - N * 3)) / (N + 7)
    P0 = -144 * (142528 + N * (23073 - N * 115))
    P1 = -12 * (6690556 + N * (2641189 + N * (245378 - N * 1495)))
    P2 = -16 * (6503401 + N * (4020794 + N * (787281 + N * (47323 - N * 385))))
    P3 =
      -(6386660 + N * (4411997 + N * (951821 + N * (61659 - N * 665)))) *
      (N + 7)
    A = P2 / (3 * P3)
    P = A * (P1 / P2 - A)
    Q = A * (2 * A * A - P1 / P3) + P0 / P3
    R = sqrt(-P^3)
    TH = acos(-Q / (2 * R)) / 3
    R = 2 * R^(1 / 3)
    TP = 2 * pi / 3
    A1 = -A + R * cos(TH)
    A2 = -A + R * cos(TH + 2 * TP)
    A3 = -A + R * cos(TH + TP)
    G[1, GMS+5] = (1 - N * A1) / NP
    G[2:NP, GMS+5] .= (1 + A1) / NP
    PTS[GMS+5] = NP
    G[1, GMS+6] = (1 - N * A2) / NP
    G[2:NP, GMS+6] .= (1 + A2) / NP
    PTS[GMS+6] = NP
    G[1, GMS+7] = (1 - N * A3) / NP
    G[2:NP, GMS+7] .= (1 + A3) / NP
    PTS[GMS+7] = NP
    W[GMS+5, IW] =
      (U7 - (A2 + A3) * U6 + A2 * A3 * U5) / (A1^2 - (A2 + A3) * A1 + A2 * A3) /
      A1^5
    W[GMS+6, IW] =
      (U7 - (A1 + A3) * U6 + A1 * A3 * U5) / (A2^2 - (A1 + A3) * A2 + A1 * A3) /
      A2^5
    W[GMS+7, IW] =
      (U7 - (A2 + A1) * U6 + A2 * A1 * U5) / (A3^2 - (A2 + A1) * A3 + A2 * A1) /
      A3^5
    G[1, GMS+8] = 4 / (N + 7)
    G[2, GMS+8] = 4 / (N + 7)
    G[3:NP, GMS+8] .= 1 / (N + 7)
    PTS[GMS+8] = NP * N / 2
    W[GMS+8, IW] = 10 * (N + 7)^6 / (729 * N6)
    G[1, GMS+9] = 11 / (N + 7) / 2
    G[2, GMS+9] = 5 / (N + 7) / 2
    G[3:NP, GMS+9] .= 1 / (N + 7)
    PTS[GMS+9] = NP * N
    W[GMS+9, IW] = 64 * (N + 7)^6 / (6561 * N6)
    W[4, IW] = W[4, IW+1]
    W[7, IW] = W[7, IW+1]
    IW = IW - 1
    G[1, 8] = 9 / (N + 9)
    G[2:NP, 8] .= 1 / (N + 9)
    PTS[8] = NP
    G[1, 9] = 7 / (N + 9)
    G[2, 9] = 3 / (N + 9)
    G[3:NP, 9] .= 1 / (N + 9)
    PTS[9] = NP * N
    G[1:2, 10] .= 5 / (N + 9)
    G[3:NP, 10] .= 1 / (N + 9)
    PTS[10] = NP * N / 2
    G[1, 11] = 5 / (N + 9)
    G[2:3, 11] .= 3 / (N + 9)
    if NP > 3
      G[4:NP, 11] .= 1 / (N + 9)
    end
    PTS[11] = NP * N * (N - 1) / 2
    W[2, IW] = -(N + 3)^9 / (6 * 256 * N6)
    W[3:4, IW] .= (N + 5)^9 / (2 * 256 * N6 * (N + 7))
    W[5:7, IW] .= -(N + 7)^9 / (256 * N8)
    W[8:11, IW] .= (N + 9)^9 / (256 * N8 * (N + 9))
    if N > 2
      G[1:4, 12] .= 3 / (N + 9)
      if NP > 4
        G[5:NP, 12] .= 1 / (N + 9)
      end
      PTS[12] = NP * N * (N - 1) * (N - 2) / 24
      W[12, IW] = W[8, IW]
    end
  end
  W[1, :] = 1 .- LinearAlgebra.transpose(PTS[2:WTS]) * W[2:WTS, :]
  NB = sum(PTS .* (W[:, 1] .* W[:, 1]))
  W[:, 2:RLS] = W[:, 2:RLS] - repeat(W[:, 1], 1, 1) * fill(1.0, 1, RLS - 1)
  W[:, 2] = W[:, 2] * sqrt(NB / sum(PTS .* W[:, 2] .* W[:, 2]))
  for K = 3:RLS
    W[:, K] =
      W[:, K] -
      W[:, 2:(K-1)] *
      LinearAlgebra.transpose(W[:, 2:(K-1)]) *
      repeat(PTS .* W[:, K], 1, 1) / NB
    W[:, K] = W[:, K] * sqrt(NB / sum(PTS .* W[:, K] .* W[:, K]))
  end
  return (G = G, W = W, PTS = PTS)
end

"""
    integratePolynomialOnSimplex(P, S)

Exact integral of a polynomial over a simplex.

# Argument
- `P`: polynomial
- `S`: simplex, given by a vector of n+1 vectors of dimension n, the simplex vertices 
"""
function integratePolynomialOnSimplex(P, S)
    gens = TypedPolynomials.variables(P)
    n = length(gens)
    if length(S) != n + 1
        error("Invalid simplex.")
    end
    for i in 1:(n+1)
        if length(S[i]) != n
            error("Invalid simplex.")
        end
    end
    v = S[n+1]    
    B = Array{Float64}(undef, n, 0)
    for i in 1:n
        B = hcat(B, S[i] - v)
    end
    Q = P(gens => v + B * vec(gens))
    s = 0.0
    for t in TypedPolynomials.terms(Q)
        coef = TypedPolynomials.coefficient(t)
        powers = TypedPolynomials.exponents(t)
        j = sum(powers)
        if j == 0
            s = s + coef
            continue
        end
        coef = coef * prod(factorial.(powers))
        s = s + coef / prod((n+1):(n+j))
    end
    return abs(LinearAlgebra.det(B)) / factorial(n) * s
end


end
