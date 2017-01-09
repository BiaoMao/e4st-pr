function [results, f, success, info, et, g, jac, xr, pimul] = e4st_solve(varargin)
%E4ST_SOLVE  Core solver for E4ST.
%
%   Examples:
%       Output argument options:
%
%       [results, f, success] = e4st_solve(...)
%       [results, f, success, info] = e4st_solve(...)
%       [results, f, success, info, et] = e4st_solve(...)
%       [results, f, success, info, et, g, jac] = e4st_solve(...)
%       [results, f, success, info, et, g, jac, xr, pimul] = e4st_solve(...)
%
%       Input argument options:
%
%       e4st_solve(mpc)
%       e4st_solve(mpc, mpopt)
%       e4st_solve(mpc, offer, contab)
%       e4st_solve(mpc, offer, contab, mpopt)
%       e4st_solve(mpc, offer, contab, A, l, u)
%       e4st_solve(mpc, offer, contab, A, l, u, mpopt)
%       e4st_solve(mpc, offer, contab, A, l, u, mpopt, N, fparm, H, Cw)
%
%   mpc is a MATPOWER case file or case struct with the fields baseMVA, bus,
%   gen, branch, and (optionally) areas. It may also include a 'contingencies'
%   field (in place of contab argument). The offer argument can be a struct
%   or a matrix. If it is a struct, it has the following fields for active
%   power quantities, each an ng x 1 vector ...
%       offer
%           .PositiveActiveReservePrice
%           .PositiveActiveReserveQuantity
%           .NegativeActiveReservePrice
%           .NegativeActiveReserveQuantity
%           .PositiveActiveDeltaPrice
%           .NegativeActiveDeltaPrice
%           .PositiveActiveReservePrice2        (optional quadratic term)
%           .NegativeActiveReservePrice2        (optional quadratic term)
%           .PositiveActiveDeltaPrice2          (optional quadratic term)
%           .NegativeActiveDeltaPrice2          (optional quadratic term)
%           .ActiveContractMin                  (optional)
%           .ActiveContractMax                  (optional)
%           .PminFactor                         (optional)
%   ... and optionally, the corresponding for reactive power ...
%           .PositiveReactiveReservePrice       (optional)
%           .PositiveReactiveReserveQuantity    (optional)
%           .NegativeReactiveReservePrice       (optional)
%           .NegativeReactiveReserveQuantity    (optional)
%           .PositiveReactiveDeltaPrice         (optional)
%           .NegativeReactiveDeltaPrice         (optional)
%           .PositiveReactiveReservePrice2      (optional quadratic term)
%           .NegativeReactiveReservePrice2      (optional quadratic term)
%           .PositiveReactiveDeltaPrice2        (optional quadratic term)
%           .NegativeReactiveDeltaPrice2        (optional quadratic term)
%           .ReactiveContractMin                (optional)
%           .ReactiveContractMax                (optional)
%   If offer is a matrix, the first ng rows contain the active power
%   quantities and the 2nd set of ng rows (optional) contain the reactive
%   power quantities. The columns correspond to the fields listed above
%   in the listed order.
%
%   Alternatively, the offer argument can be omitted and the fields
%   'reserve', 'energy_delta_cost' and 'contract' included in mpc.
%   In this case, the 'reserve', 'energy_delta_cost' and 'contract' fields
%   take the following form, where offerp refers to the first ng rows of
%   the corresponding offer matrix and offerq to the optional 2nd set of
%   ng rows:
%       .reserve
%           .cost
%               .Rp_pos     [ offerp(:, 1) ]
%               .Rp_neg     [ offerp(:, 3) ]
%               .Rp_pos2    [ offerp(:, 7) ]    (optional quadratic term)
%               .Rp_neg2    [ offerp(:, 8) ]    (optional quadratic term)
%               .Rq_pos     [ offerq(:, 1) ]    (optional)
%               .Rq_neg     [ offerq(:, 3) ]    (optional)
%               .Rq_pos2    [ offerq(:, 7) ]    (optional quadratic term)
%               .Rq_neg2    [ offerq(:, 8) ]    (optional quadratic term)
%           .cap
%               .Rp_pos     [ offerp(:, 2) ]
%               .Rp_neg     [ offerp(:, 4) ]
%               .Rq_pos     [ offerq(:, 2) ]    (optional)
%               .Rq_neg     [ offerq(:, 4) ]    (optional)
%       .energy_delta_cost
%           .dP_pos         [ offerp(:, 5) ]
%           .dP_neg         [ offerp(:, 6) ]
%           .dP_pos2        [ offerp(:, 9) ]    (optional quadratic term)
%           .dP_neg2        [ offerp(:, 10)]    (optional quadratic term)
%           .dQ_pos         [ offerq(:, 5) ]    (optional)
%           .dQ_neg         [ offerq(:, 6) ]    (optional)
%           .dQ_pos2        [ offerq(:, 9) ]    (optional quadratic term)
%           .dQ_neg2        [ offerq(:, 10)]    (optional quadratic term)
%       .contract                               (optional)
%           .Pc_min         [ offerp(:, 11)]    (optional)
%           .Pc_max         [ offerp(:, 12)]    (optional)
%           .Qc_min         [ offerq(:, 11)]    (optional)
%           .Qc_max         [ offerq(:, 12)]    (optional)
%       .pmin_factor        [ offerp(:, 13)]    (optional)
%
%   An optional 'availability_factor' field can be used to specify an
%   availability factor for each generator in each scenario. It can be
%   either an ng x 1 vector or ng x (nc+1) matrix, where each element is
%   between 0 and 1.
%
%   contab is the contingency table, type 'help apply_changes' and
%   'help idx_ct' for details about the format.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.

%   E4ST
%   Copyright (c) 2000-2017 by Power System Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.

% SECTION 1: ARGUMENT PARSING, REORDERING, UNPACKING

[baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
    Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q, ...
    caplim_map, caplim_max, caplim_min, avail_fac, ...
    toc_map, toc_cap, toc_coeff, toc_type] = e4st_args(varargin{:});

if size(N, 1) > 0
  if size(N, 1) ~= size(fparm, 1) || size(N, 1) ~= size(H,1) || ...
     size(N, 1) ~= length(Cw)
    error('e4st_solve.m: wrong internal dimensions in generalized cost parameters');
  end
end

if nargout > 5
    mpopt = mpoption(mpopt, 'opf.return_raw_der', 1);
end

% options
OUT_ALL = mpopt.out.all;
dc      = strcmp(upper(mpopt.model), 'DC');
if isfield(mpopt, 'sopf') && isfield(mpopt.sopf, 'force_Pc_eq_P0')
    FORCE_PC_EQ_P0 = mpopt.sopf.force_Pc_eq_P0;
else
    FORCE_PC_EQ_P0 = 0;     %% off by default
end

if mpopt.verbose > 0
  v = e4st_ver('all');
  fprintf('\nE4ST Version %s, %s\n', v.Version, v.Date);
  fprintf('Engineering, Economic, and Environmental Electricity Simulation Tool\n');
end

if dc                   %% force HAVE_Q to false for DC runs
    HAVE_Q = 0;
end

% Load column indices for case tables.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;
c = idx_dcline;

% If tables do not have multiplier/extra columns, append zero cols.
% Update whenever the data format changes!
if size(bus,2) < MU_VMIN
  bus = [bus zeros(size(bus,1),MU_VMIN-size(bus,2)) ];
end
if size(gen,2) < MU_QMIN
  gen = [ gen zeros(size(gen,1),MU_QMIN-size(gen,2)) ];
end
if size(branch,2) < MU_ANGMAX
  branch = [ branch zeros(size(branch,1),MU_ANGMAX-size(branch,2)) ];
end
if ~isempty(dcline) && size(dcline,2) < c.MU_QMAXT
  dcline = [ dcline zeros(size(dcline,1),c.MU_QMAXT-size(dcline,2)) ];
end

% get number of contingencies & list of labels
[xx, ii] = sort(contab(:, CT_LABEL)); %sort in ascending contingency label
contab = contab(ii, :);
clist = unique(contab(:, CT_LABEL));

% Filter-out uncomitted equipment from data; the row indices in
% contab must be modified accordingly.  This might lead to
% deleting contingencies... yuk... must track this too if this
% is to be helpful to a higher-level decommitment routine...

% More on higher-level unit de-commitment: an outer de-commit
% decision for a generator or a line might make a corresponding
% contingency in the original list irrelevant; the problem becomes
% entirely different, with different base case probability and so on.
% To keep things consistent, on output, the post-contingency
% flows for deleted contingencies should still have a numbered spot
% in the data but be empty, which must be correctly interpreted
% by the calling unit de-commitment code.

base_on_gen  = find(gen(:, GEN_STATUS) >  0);
base_off_gen = find(gen(:, GEN_STATUS) <= 0);
base_on_branch  = find(branch(:, BR_STATUS) >  0);
base_off_branch = find(branch(:, BR_STATUS) <= 0);
gen_original = gen;                 % save these three
branch_original = branch;
clist_original = clist;
gen =  gen(base_on_gen, :);         % now these contain only equipment
branch = branch(base_on_branch, :); % committed on input
ng = size(gen_original, 1);         % original size(gen), at least temporarily
if size(gencost,1) == ng            % while we do something with it
  gencost = gencost(base_on_gen, :);
else
  gencost = gencost( [base_on_gen; base_on_gen+ng], :);
end
if ~isempty(dcline)
  base_on_dcline  = find(dcline(:, c.BR_STATUS) >  0);
  base_off_dcline = find(dcline(:, c.BR_STATUS) <= 0);
  dcline_original = dcline;
  dcline = dcline(base_on_dcline, :);
end
if ~isempty(iflims)
  ifmap = iflims.map;
  ifmap_original = ifmap;
  % update branch indices in ifmap, remove lines that are out-of-service
  e2i = zeros(size(branch_original, 1), 1);
  e2i(base_on_branch) = (1:size(branch, 1))';
  d = sign(ifmap(:, 2));
  br = abs(ifmap(:, 2));
  ifmap(:, 2) = d .* e2i(br);
  ifmap(ifmap(:, 2) == 0, :) = [];  %% delete branches that are out
end
if ~isempty(caplim_map)
  caplim_map = caplim_map(:, base_on_gen);
end
if ~isempty(avail_fac)
  avail_fac = avail_fac(base_on_gen, :);        % reduce to committed gens
  af_nc = size(avail_fac, 2);
  if af_nc > 1 && af_nc ~= length(clist) + 1    % check for proper # of cols
    error('e4st_solve.m: # of columns in mpc.availability_factor (%d) must equal 1 or # contingencies + 1 (%d)', ...
        af_nc, length(clist) + 1);
  end
end
if ~isempty(toc_map)
  toc_map = toc_map(:, base_on_gen, :);
  toc_coeff = toc_coeff(base_on_gen, :);
end

% Do the same to offer table/struct
% offer_original = offer;
if HAVE_Q
  offer = offer([base_on_gen; base_on_gen+ng], :);
else
  offer = offer(base_on_gen, :);
end

% Look for contingency table rows that act on branches or generators that
% are turned off on input (and therefore have now been removed). If so, delete
% the corresponding contingency table row. At the end, see if a contingency
% label in clist points to (now) nonexistent labels in the contingency
% table and if so entirely delete the contingency label from clist.
% NOTE: Does not catch changes that act on all rows using a row value of 0,
% or e.g. all gens in an area... yet ***************

rowcomlist = ones(size(contab,1), 1);
for l = 1:size(contab, 1)
  if contab(l, CT_ROW) ~= 0 && ...
        ((contab(l, CT_TABLE) == CT_TGEN && ...
            gen_original(contab(l, CT_ROW), GEN_STATUS) <= 0) || ...
         (contab(l, CT_TABLE) == CT_TBRCH && ...
            branch_original(contab(l, CT_ROW), BR_STATUS) <= 0))
    rowcomlist(l) = 0;
  end
end
contab = contab(rowcomlist ~= 0, :);
clabelcomlist = ones(size(clist));
for l = 1:length(clist)
  if ~any(clist(l) == contab(:, CT_LABEL))
    clabelcomlist(l) = 0;
  end
end
clist = clist(clabelcomlist ~= 0);
% clabeldecomlist = clist(clabelcomlist == 0); % contingencies with these labels were thrown out
nc = length(clist); % redefined

% resize avail_fac to make it ng x (nc+1)
if ~isempty(avail_fac)
  if size(avail_fac, 2) > 1
    avail_fac = avail_fac(:, [1; clabelcomlist] ~= 0);  % reduce for removed contingencies
  else
    avail_fac = avail_fac * ones(1, nc+1);  % expand vector to matrix
  end
else
  avail_fac = ones(ng, nc+1);               % all ones if not given
end

% Now catch renumbering of CT_ROW in contab due to deletion of rows
% for equipment that came with STATUS=off on input. Remember that at this
% point contab points only to rows that denote equipment that is to be kept
% in the analysis.

newrowgen = cumsum(gen_original(:,GEN_STATUS) > 0);   % newrowgen(4) contains
newrowbrch = cumsum(branch_original(:,BR_STATUS) > 0);% new row index for gen 4,
if size(gencost,1) == ng                              % assuming gen4 is kept
  newrowgencost = newrowgen;
else
  newrowgencost = cumsum([  gen_original(:,GEN_STATUS);
                            gen_original(:,GEN_STATUS) ]  > 0);
end
for l = 1:size(contab,1)
  if contab(l, CT_TABLE) == CT_TGEN && contab(l, CT_ROW) > 0
    contab(l, CT_ROW) = newrowgen(contab(l, CT_ROW));
  elseif contab(l, CT_TABLE) == CT_TGENCOST && contab(l, CT_ROW) > 0
    contab(l, CT_ROW) = newrowgencost(contab(l, CT_ROW));
  elseif contab(l, CT_TABLE) == CT_TBRCH && contab(l, CT_ROW) > 0
    contab(l, CT_ROW) = newrowbrch(contab(l, CT_ROW));
  end
end

% compute probability of nominal (base) flow as the complement of the
% sum of the probabilities of the contingencies
prob0 = 1;
for i=1:nc
  tmp = contab(clist(i) == contab(:, CT_LABEL), CT_PROB);
  prob0 = prob0 - tmp(1);
end

% Renumber buses consecutively (rows themselves in bus() are not permuted)
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
if ~isempty(dcline)
  e2i = sparse(max(i2e), 1);
  e2i(i2e) = (1:size(bus, 1))';
  dcline(:, c.F_BUS) = e2i( dcline(:, c.F_BUS) );
  dcline(:, c.T_BUS) = e2i( dcline(:, c.T_BUS) );
end

% Sort generators in order of increasing bus number; rows _are_ permuted.
% Therefore any contingency including a change to a gen() table row
% must be adjusted accordingly.  Offers and gencosts must also be permuted.
ng = size(gen,1);
[tmp, igen] = sort(gen(:, GEN_BUS));
[tmp, inv_gen_ord] = sort(igen);  % save for inverse reordering at the end
gen  = gen(igen, :);
if ng == size(gencost,1)
  gencost = gencost(igen, :);
  inv_gencost_ord = inv_gen_ord;
else
  gencost = gencost( [igen; igen+ng], :);
  [tmp, inv_gencost_ord] = sort([igen; igen+ng]);
end
[gencost_p, gencost_q] = pqcost(gencost, ng);
if ~isempty(caplim_map)
  caplim_map = caplim_map(:, igen);
end
avail_fac = avail_fac(igen, :);
if ~isempty(toc_map)
  toc_map = toc_map(:, igen, :);
  toc_coeff = toc_coeff(igen, :);
end
igenmods = find(contab(:,CT_TABLE) == CT_TGEN & contab(:,CT_ROW) ~= 0);%which rows modify gen table
if ~isempty(igenmods)
  contab(igenmods, CT_ROW) = inv_gen_ord(contab(igenmods, CT_ROW));
end
igencostmods = find(contab(:,CT_TABLE) == CT_TGENCOST & contab(:,CT_ROW) ~= 0);%which rows modify gencost table
if ~isempty(igencostmods)
  contab(igencostmods, CT_ROW) = inv_gencost_ord(contab(igencostmods, CT_ROW));
end
% Offer processing - from here on we have unified source of offer data
PositiveActiveReservePrice      = offer(igen, 1);
PositiveActiveReserveQuantity   = offer(igen, 2);
NegativeActiveReservePrice      = offer(igen, 3);
NegativeActiveReserveQuantity   = offer(igen, 4);
PositiveActiveDeltaPrice        = offer(igen, 5);
NegativeActiveDeltaPrice        = offer(igen, 6);
PositiveActiveReservePrice2     = offer(igen, 7);
NegativeActiveReservePrice2     = offer(igen, 8);
PositiveActiveDeltaPrice2       = offer(igen, 9);
NegativeActiveDeltaPrice2       = offer(igen, 10);
ActiveContractMin               = offer(igen, 11);
ActiveContractMax               = offer(igen, 12);
if size(offer, 2) >= 13
    PminFactor                  = offer(igen, 13);
else
    PminFactor                  = zeros(length(igen), 1);
end
if HAVE_Q
  PositiveReactiveReservePrice      = offer(ng+igen, 1);
  PositiveReactiveReserveQuantity   = offer(ng+igen, 2);
  NegativeReactiveReservePrice      = offer(ng+igen, 3);
  NegativeReactiveReserveQuantity   = offer(ng+igen, 4);
  PositiveReactiveDeltaPrice        = offer(ng+igen, 5);
  NegativeReactiveDeltaPrice        = offer(ng+igen, 6);
  PositiveReactiveReservePrice2     = offer(ng+igen, 7);
  NegativeReactiveReservePrice2     = offer(ng+igen, 8);
  PositiveReactiveDeltaPrice2       = offer(ng+igen, 9);
  NegativeReactiveDeltaPrice2       = offer(ng+igen, 10);
  ReactiveContractMin               = offer(ng+igen, 11);
  ReactiveContractMax               = offer(ng+igen, 12);
end
if any(any(offer(:, 7:10))) % we have quadratic inc/dec/reserve cost terms
  HAVE_QUADRATIC = 1;
else
  HAVE_QUADRATIC = 0;
end
if FORCE_PC_EQ_P0 && any(any(isfinite(offer(:, 11:12))))
  fprintf('\ne4st_solve: WARNING: Using limits on Pc (or Qc) is not recommended\n');
  fprintf('                 when the ''sopf.force_Pc_eq_P0'' option is enabled.\n\n');
end

% SECTION 2: CONTINGENCY PROCESSING AND SYSTEM AUGMENTATION

% The big-system building must be a two-pass process since
% start/end indices of both variables and constraints cannot be known
% until all contingencies are processed and we see how many
% injections and branches are active in each flow.  After the first
% pass, indices can be computed and then additional linear constraints
% can be built in the second pass.

% FIRST PASS. Augmented network information is built, as well as
% basic augmented generator and cost tables.
ng(1) = size(gen, 1);
nb(1) = size(bus, 1);   % not supposed to change, but..
nl(1) = size(branch, 1);
ndc(1) = size(dcline, 1);
if ~isempty(iflims)
    nif(1) = size(iflims.lims, 1);
    nifm(1) = size(ifmap, 1);
    ifidmax(1) = max(ifmap(:, 1));
else
    nif(1) = 0;
    nifm(1) = 0;
    ifidmax(1) = 0;
end
Augmented_bus = bus;
Augmented_branch = branch;
Augmented_gen = gen;
Augmented_gencost_p = modcost(gencost_p, prob0);
if ~isempty(gencost_q)
  Augmented_gencost_q = modcost(gencost_q, prob0);
else
  Augmented_gencost_q = [];
end
if ~isempty(dcline)
  Augmented_dcline = dcline;
end
if ~isempty(iflims)
    Augmented_ifmap = ifmap;
    Augmented_iflims = iflims.lims;
end
gen_stat = [];
branch_stat = [];
% dcline_stat = [];
for k = 1:nc
  clabel = clist(k);
  kthcontab = contab(clabel == contab(:, CT_LABEL), :); % get rows addressing kth contingency
  ng(k+1) = ng(1); % initialize to 'full size' and we'll decrement from there
  nb(k+1) = nb(1);
  nl(k+1) = nl(1);
  ndc(k+1) = ndc(1);
  nif(k+1) = nif(1);
  nifm(k+1) = nifm(1);
  ifidmax(k+1) = ifidmax(k) + ifidmax(1);

  % apply the modifications
  mpck = struct('bus', bus, 'gen', gen, 'branch', branch, 'gencost', gencost);
  mpck = apply_changes(clabel, mpck, contab);

  % working copies of tables for k-th contingency
  [buswork, genwork, branchwork] = deal(mpck.bus, mpck.gen, mpck.branch);
  [gencost_p_work, gencost_q_work] = pqcost(mpck.gencost, ng(1));

  % now catch any dimension changes that result from them and construct
  % augmented network data for this contingency
  nb(k+1) = size(bus,1);
  busbastmp = sum(nb(1:k));  % total # buses in base and previous flows
  buswork(:, BUS_I) = bus(:, BUS_I) + busbastmp;
  % augmented problem has multiple REF buses
  % opf() has been updated to fix all such REF bus voltage angles and
  % emit a warning in case that was not intentional
  gen_stat = [ gen_stat genwork(:, GEN_STATUS) ];
  ii = find(genwork(:, GEN_STATUS) <= 0);
  genwork(ii, : ) = [];
  gencost_p_work(ii, :) = [];
  if ~isempty(gencost_q_work)
    gencost_q_work(ii, : ) = [];
  end
  ng(k+1) = size(genwork, 1);
  genwork(:, GEN_BUS) = genwork(:, GEN_BUS) + busbastmp; % bus #s kth island
  branch_stat = [ branch_stat branchwork(:,BR_STATUS) ];
%  branchwork(:, RATE_A) = branchwork(:, RATE_B);
  branchwork(branchwork(:, BR_STATUS) <= 0, : ) = [];
  nl(k+1) = size(branchwork, 1);
  brbastmp = sum(nl(1:k));  % total # branches in base and previous flows
  branchwork(:, F_BUS) = branchwork(:, F_BUS) + busbastmp;
  branchwork(:, T_BUS) = branchwork(:, T_BUS) + busbastmp;
  if ~isempty(dcline)
    dclinework = dcline;    %% assume no changes per contingency (for now)
%     dcline_stat = [ dcline_stat dclinework(:,BR_STATUS) ];
%     dclinework(dclinework(:, BR_STATUS) <= 0, : ) = [];
%     ndc(k+1) = size(dclinework, 1);
    dclinework(:, c.F_BUS) = dclinework(:, c.F_BUS) + busbastmp;
    dclinework(:, c.T_BUS) = dclinework(:, c.T_BUS) + busbastmp;
    Augmented_dcline = [Augmented_dcline; dclinework];
  end
  if ~isempty(iflims)
    ifmapwork = ifmap;          %% assume no changes per contingency
    iflimswork = iflims.lims;   %% assume no changes per contingency
    ifbastmp = ifidmax(k);      %% max if id num in previous flows
    ifmapwork(:, 1)  = ifmapwork(:, 1) + ifbastmp;
    iflimswork(:, 1) = iflimswork(:, 1) + ifbastmp;
    e2i = zeros(size(mpck.branch, 1), 1);
    e2i(mpck.branch(:, BR_STATUS) >  0) = (1:size(branchwork, 1))';
    d = sign(ifmapwork(:, 2));
    br = abs(ifmapwork(:, 2));
    ifmapwork(:, 2) = d .* (e2i(br) + brbastmp);
    ifmapwork(ifmapwork(:, 2) == brbastmp, :) = [];  %% delete branches that are out
    nifm(k+1) = size(ifmapwork, 1);
    Augmented_ifmap = [Augmented_ifmap; ifmapwork];
    Augmented_iflims = [Augmented_iflims; iflimswork];
  end
  Augmented_bus = [Augmented_bus; buswork];
  Augmented_gen = [Augmented_gen; genwork];
  Augmented_branch = [Augmented_branch; branchwork];
  Augmented_gencost_p = [Augmented_gencost_p; ...
                         modcost(gencost_p_work, kthcontab(1, CT_PROB)) ];
  Augmented_gencost_q = [Augmented_gencost_q; ...
                         modcost(gencost_q_work, kthcontab(1, CT_PROB)) ];
end % for k=contingencies 

Augmented_gencost = [ Augmented_gencost_p; Augmented_gencost_q];

% SECTION 3: INDEXING SCHEME FOR VARIABLES AND NONLINEAR CONSTRAINTS

% build variable indexing scheme
nb_total = sum(nb(:));
ng_total = sum(ng(:));
nl_total = sum(nl(:));
% ndc_total = sum(ndc(:));
thbas(1) = 1;                      thend(1) = thbas(1) + nb(1) - 1;
vbas(1)  = thbas(1) + nb_total;    vend(1)  = vbas(1)  + nb(1) - 1;
pgbas(1) = 2*nb_total + 1;         pgend(1) = pgbas(1) + ng(1) - 1;
qgbas(1) = 2*nb_total+ng_total+1;  qgend(1) = qgbas(1) + ng(1) - 1;
if nc > 0
  for k = 2:nc+1
    thbas(k) = thend(k-1) + 1;     thend(k) = thbas(k) + nb(k) - 1;
    vbas(k)  = vend(k-1)  + 1;     vend(k)  = vbas(k)  + nb(k) - 1;
    pgbas(k) = pgend(k-1) + 1;     pgend(k) = pgbas(k) + ng(k) - 1;
    qgbas(k) = qgend(k-1) + 1;     qgend(k) = qgbas(k) + ng(k) - 1;
  end
end
% contract active energy quantities
pcbas = qgend(nc+1) + 1;
pcend = pcbas + ng(1) - 1;
% upward and downward P reserves
rPpbas = pcend + 1;                rPpend = rPpbas + ng(1) - 1; % positive P reserve
rPmbas = rPpbas + ng(1);           rPmend = rPmbas + ng(1) - 1; % downward P reserve
% P deviations from contracted base injections (nc+1 sets of them, one for each flow)
dPpbas(1) = rPmend + 1;             dPpend(1) = dPpbas(1) + ng(1) - 1;
dPmbas(1) = dPpbas(1)+ng_total;     dPmend(1) = dPmbas(1) + ng(1) - 1;
if nc > 0
  for k = 2:nc+1
    dPpbas(k) = dPpend(k-1) + 1;    dPpend(k) = dPpbas(k) + ng(k) - 1;
    dPmbas(k) = dPmend(k-1) + 1;    dPmend(k) = dPmbas(k) + ng(k) - 1;
  end
end

% if reactive stuff has prices, generate pointers to reactive contract
% quantities, reserves and deviations in the vector of optimization variables
if HAVE_Q
  qcbas = dPmend(nc) + 1;        qcend = qcbas + ng(1) - 1;
  rQpbas = qcend + 1;            rQpend = rQpbas + ng(1) - 1; % upward Q reserve
  rQmbas = rQpend + 1;           rQmend = rQmbas + ng(1) - 1; % downward Q reserve
  % Q deviations from base flow
  dQpbas(1) = rQmend + 1;        dQpend(1) = dQpbas(1) + ng(1) - 1;
  dQmbas(1) = dQpbas(1)+ng_total; dQmend(1) = dQmbas(1) + ng(1) - 1;
  if nc > 0
    for k = 2:nc+1
      dQpbas(k) = dQpend(k-1) + 1;   dQpend(k) = dQpbas(k) + ng(k) - 1;
      dQmbas(k) = dQmend(k-1) + 1;   dQmend(k) = dQmbas(k) + ng(k) - 1;
    end
  end
  nvars = dQmend(nc+1);  % last variable when reactive portion is considered
else
  nvars = dPmend(nc+1);  % last variable when reactive portion not considered.
end

% Reorder the user-provided matrix N, which assumes a variable ordering
% (nc+1)*nb angles, 
% (nc+1)*nb voltage magnitudes, 
% (nc+1)*ng activeinjections (including data for those not committed),
% (nc+1)*ng reactive injections, 
% ng active contract quantities, 
% ng R+ vars,
% ng R- vars, 
% (nc+1)*ng dP+ vars, 
% (nc+1)*ng dP- vars, 
% ng reactive contract quantities, 
% ng Rq+ vars, 
% ng Rq- vars, 
% (nc+1)*ng dQ+ vars,
% (nc+1)*ng dQ- vars
% for a total of 2*(nc+1)*nb + 6*(nc+1)*ng +6*ng variables

ng_orig = size(gen_original, 1);
if HAVE_Q
  nvars_orig = 2*(nc+1)*nb(1) + 6*(nc+1)*ng_orig +6*ng_orig;
else
  nvars_orig = 2*(nc+1)*nb(1) + 4*(nc+1)*ng_orig +3*ng_orig;
end
if ~isempty(N)
  if size(N,2) < nvars_orig
    if issparse(N)
      N(:, size(N,2)+1:nvars_orig) = sparse(size(N,1),nvars_orig-size(N,2));
    else
      N(:, size(N,2)+1:nvars_orig) = zeros(size(N,1),nvars_orig-size(N,2));
    end
  end
  % Now delete locations that refer to units that were not committed in
  % the original data; to do this, construct a vector of indices to
  % those locations and use it to assign an empty matrix to the
  % corresponding columns of N.
  l = pgbas(1);  % this is still the location of the first injection
  del_locs = l-1+base_off_gen;                   % pi0
  l = l + ng_orig;
  if nc > 1
    for k = 2:nc+1
      del_locs = [ del_locs; l-1+base_off_gen ]; % pik
      l = l + ng_orig;
    end
  end
  del_locs = [ del_locs; l-1+base_off_gen ];     % qi0
  l = l + ng_orig;
  if nc > 1
    for k = 2:nc+1
      del_locs = [ del_locs; l-1+base_off_gen ]; % qik
      l = l + ng_orig;
    end
  end
  del_locs = [ del_locs; l-1+base_off_gen ];     % pc
  l = l + ng_orig;
  del_locs = [ del_locs; l-1+base_off_gen ];     % rPp
  l = l + ng_orig;
  del_locs = [ del_locs; l-1+base_off_gen ];     % rPm
  l = l + ng_orig;
  del_locs = [ del_locs; l-1+base_off_gen ];     % dPp0
  l = l + ng_orig;
  if nc > 1
    for k = 2:nc+1
      del_locs = [ del_locs; l-1+base_off_gen ]; % dPpik
      l = l + ng_orig;
    end
  end
  del_locs = [ del_locs; l-1+base_off_gen ];     % dPm0
  l = l + ng_orig;  
  if nc > 1
    for k = 2:nc+1
      del_locs = [ del_locs; l-1+base_off_gen ]; % dPmik
      l = l + ng_orig;
    end
  end
  if HAVE_Q
    del_locs = [ del_locs; l-1+base_off_gen ];   % qc
    l = l + ng_orig;
    del_locs = [ del_locs; l-1+base_off_gen ];   % rQp
    l = l + ng_orig;
    del_locs = [ del_locs; l-1+base_off_gen ];   % rQm
    l = l + ng_orig;
    del_locs = [ del_locs; l-1+base_off_gen ];   % dQp0
    l = l + ng_orig;
    if nc > 1
      for k = 2:nc+1
        del_locs = [ del_locs; l-1+base_off_gen ]; % dQpik
        l = l + ng_orig;
      end
    end
    del_locs = [ del_locs; l-1+base_off_gen ];   % dQm0
    l = l + ng_orig;
    if nc > 1
      for k = 2:nc+1
        del_locs = [ del_locs; l-1+base_off_gen ]; % dQmik
        l = l + ng_orig;
      end
    end
  end
  % Ready; annihilate columns
  N(:, del_locs) = [];
  % Now reorder sets of columns according to internal generator order
  % (we still have columns related to nonexistent variables, such as
  % the injection of an ousted generator in a contingency); those cols
  % are to be annihilated later).
  l = pgbas(1);  % this is still the location of the first injection
  N(:, l:l+ng(1)-1) =  N(:, l-1+igen);           % pi0
  l = l + ng(1);
  if nc > 1
    for k = 2:nc+1
      N(:, l:l+ng(1)-1) =  N(:, l-1+igen);       % pik
      l = l + ng(1);
    end
  end
  N(:, l:l+ng(1)-1) =  N(:, l-1+igen);           % qi0
  l = l + ng(1);
  if nc > 1
    for k = 2:nc+1
      N(:, l:l+ng(1)-1) =  N(:, l-1+igen);       % qik
      l = l + ng(1);
    end
  end
  N(:, l:l+ng(1)-1) =  N(:, l-1+igen);           % pc
  l = l + ng(1);
  N(:, l:l+ng(1)-1) =  N(:, l-1+igen);           % rPp
  l = l + ng(1);
  N(:, l:l+ng(1)-1) =  N(:, l-1+igen);           % rPm
  l = l + ng(1);
  N(:, l:l+ng(1)-1) =  N(:, l-1+igen);           % dPp0
  l = l + ng(1);
  if nc > 1
    for k = 2:nc+1
      N(:, l:l+ng(1)-1) =  N(:, l-1+igen);       % dPpik
      l = l + ng(1);
    end
  end
  N(:, l:l+ng(1)-1) =  N(:, l-1+igen);           % dPm0
  l = l + ng(1);  
  if nc > 1
    for k = 2:nc+1
      N(:, l:l+ng(1)-1) =  N(:, l-1+igen);       % dPmik
      l = l + ng(1);
    end
  end
  if HAVE_Q
    N(:, l:l+ng(1)-1) =  N(:, l-1+igen);         % qc
    l = l + ng(1);
    N(:, l:l+ng(1)-1) =  N(:, l-1+igen);         % rQp
    l = l + ng(1);
    N(:, l:l+ng(1)-1) =  N(:, l-1+igen);         % rQm
    l = l + ng(1);
    N(:, l:l+ng(1)-1) =  N(:, l-1+igen);         % dQp0
    l = l + ng(1);
    if nc > 1
      for k = 2:nc+1
        N(:, l:l+ng(1)-1) =  N(:, l-1+igen);     % dQpik
        l = l + ng(1);
      end
    end
    N(:, l:l+ng(1)-1) =  N(:, l-1+igen);         % dQm0
    l = l + ng(1);
    if nc > 1
      for k = 2:nc+1
        N(:, l:l+ng(1)-1) =  N(:, l-1+igen);     % dQmik
        l = l + ng(1);
      end
    end
  end
  % Now all that remains is to eliminate columns for generators
  % ousted in contingencies.  Again, build vector of indices of
  % such columns.
  del_locs = [];
  l = pgbas(2);  % should point to base of pik for contingencies
  if nc > 1
    for k = 2:nc
      del_locs = [ del_locs; l-1+find(gen_stat(:,k-1) <= 0)];
      l = l + ng(1);
    end
  end
  l = l + ng(1); % now points to qgbas(2) in the order in which N is so far
  if nc > 1
    for k = 2:nc
      del_locs = [ del_locs; l-1+find(gen_stat(:,k-1) <= 0)];
      l = l + ng(1);
    end
  end
  l = l + 4*ng(1);  % Skip pc, rPp, rPm, dPp0; now it points to dPp1
  if nc > 1
    for k = 2:nc
      del_locs = [ del_locs; l-1+find(gen_stat(:,k-1) <= 0)];
      l = l + ng(1);
    end
  end
  l = l + ng(1);    % Skip dPm0
  if nc > 1
    for k = 2:nc
      del_locs = [ del_locs; l-1+find(gen_stat(:,k-1) <= 0)];
      l = l + ng(1);
    end
  end
  if HAVE_Q
    l = l + 4*ng(1); % Skip qc, rQp, rQm, dQp0; now it points to dQp1
    if nc > 1
      for k = 2:nc
        del_locs = [ del_locs; l-1+find(gen_stat(:,k-1) <= 0)];
        l = l + ng(1);
      end
    end
    l = l + ng(1);   % Skip dQm0
    if nc > 1
      for k = 2:nc
        del_locs = [ del_locs; l-1+find(gen_stat(:,k-1) <= 0)];
        l = l + ng(1);
      end
    end
  end
  N(:, del_locs) = [];
end % if ~isempty(N)

% build constraint indexing scheme for nonlinear constraints
pmsmbas(1) = 1;                    pmsmend(1) = pmsmbas(1) + nb(1) - 1;
qmsmbas(1) = pmsmbas(1)+ nb_total; qmsmend(1) = qmsmbas(1) + nb(1) - 1;
sfbas(1) = 2*nb_total + 1;         sfend(1) = sfbas(1) + nl(1) - 1;
stbas(1) = sfbas(1) + nl_total;    stend(1) = stbas(1) + nl(1) - 1;
if nc > 0
  for k = 2:nc+1
    pmsmbas(k) = pmsmend(k-1) + 1; pmsmend(k) = pmsmbas(k) + nb(k) - 1;
    qmsmbas(k) = qmsmend(k-1) + 1; qmsmend(k) = qmsmbas(k) + nb(k) - 1;
    sfbas(k) = sfend(k-1) + 1;     sfend(k) = sfbas(k) + nl(k) - 1;
    stbas(k) = stend(k-1) + 1;     stend(k) = stbas(k) + nl(k) - 1;
  end
end

% SECTION 4: LINEAR CONSTRAINTS FOR ACTIVE POWER RESERVES AND INCREMENTS

aug_gen_stat = [ ones(ng(1),1) gen_stat]; % table includes base case column
                                          % as opposed to gen_stat

% Start defining the constraints.  First define ng(1) upward P reserve variables
% rpp <= rppmax, but do it backwards so we get multipliers with
% correct sign; MINOS writes a Lagrangian F(X) - lambda^T * G(X),
% so we want the lower limit to be usually binding to get positive
% multipliers, so we rewrite this as -rppmax <= -rpp
A1 = sparse((1:ng(1))', (rPpbas:rPpend)', -ones(ng(1),1), ng(1), nvars);
%l1 = -min(PositiveActiveReserveQuantity, gen(:, RAMP_10)) / baseMVA;
l1 = -PositiveActiveReserveQuantity / baseMVA;
u1 = Inf(ng(1),1);
lc1bas = stend(nc+1)+size(Au,1) + 1;    % linear constraint set 1 start index
lc1end = lc1bas + ng(1)-1;              % linear constraint set 1 end index
% Now define ng(1) downward P reserve variables via rpm <= rpmmax,
% or, for the reasons explained above, as -rpmmax <= -rpm
A2 = sparse((1:ng(1))', (rPmbas:rPmend)', -ones(ng(1),1), ng(1), nvars);
%l2 = -min(NegativeActiveReserveQuantity, gen(:, RAMP_10)) / baseMVA;
l2 = -NegativeActiveReserveQuantity / baseMVA;
u2 = Inf(ng(1), 1);
lc2bas = lc1end + 1;
lc2end = lc2bas + ng(1) - 1;
% The actual deviations from base flow must not exceed physical ramp rates
% we'll get negative multiplier for right bound, fix when picking up lambdas
A3 = sparse(0,0); l3 = []; u3 = [];
for k = 2:nc+1
  ramp10k = Augmented_gen(sum(ng(1:k-1)) + (1:ng(k)), RAMP_10);
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A3 = [ A3 ;
         sparse([ 1:ng(k)           1:ng(k)]', ...
                [ pgbas(1)-1+ii;    (pgbas(k):pgend(k))' ], ...
                [ -ones(ng(k),1);   ones(ng(k),1)        ], ...
                ng(k), nvars)
        ];
  l3 = [ l3;
         -ramp10k/baseMVA ];
  u3 = [ u3;
         ramp10k/baseMVA  ];
end
lc3bas = lc2end + 1;
lc3end = lc3bas + size(A3,1) - 1;
% alpha-controlled equality of Pi0 and Pc
if FORCE_PC_EQ_P0
  alpha = 0;
else
  alpha = 2*max(max(abs(gen(:,[PMIN PMAX QMIN QMAX]))))/baseMVA+10;
end
A4 = sparse([1:ng(1) 1:ng(1)]', [pgbas(1):pgend(1) pcbas:pcend]', ...
            [ones(ng(1),1); -ones(ng(1),1)],   ng(1), nvars);
l4 = -alpha*ones(ng(1),1);
u4 = alpha*ones(ng(1),1);
lc4bas = lc3end + 1;
lc4end = lc4bas + ng(1) - 1;
% bounds on Pc
A4B = sparse([1:ng(1)]', [pcbas:pcend]', ones(ng(1),1), ng(1), nvars);
l4B = ActiveContractMin / baseMVA;
u4B = ActiveContractMax / baseMVA;
lc4Bbas = lc4end + 1;
lc4Bend = lc4Bbas + ng(1) - 1;
% capacity limits
if ~isempty(caplim_map)
  nclb = size(caplim_map, 1);     %% number of caplim bounds
  A4C = sparse(nclb, nvars);
  A4C(:, (rPpbas:rPpend)') = caplim_map;
%   [ii, jj, ss] = find(caplim_map);
%   A4C = sparse(ii, jj + rPpbas-1, ss, nclb, nvars);
  if isempty(caplim_min)
    l4C = -Inf(nclb, 1);
  else
    l4C = caplim_min / baseMVA;
  end
  if isempty(caplim_max)
    u4C = Inf(nclb, 1);
  else
    u4C = caplim_max / baseMVA;
  end
else
  nclb = 0;
  A4C = sparse(nclb, nvars);
  l4C = [];
  u4C = [];
end
lc4Cbas = lc4Bend + 1;
lc4Cend = lc4Cbas + nclb - 1;

% The increment variables from the base case for generators active in kth
% flow cannot be negative. 0 <= dPp
A5 = sparse(0,0); l5 = []; u5 = [];
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A5 = [ A5;
        sparse((1:ng(k))', (dPpbas(k):dPpend(k))', ones(ng(k),1), ng(k), nvars) ];
  l5 = [ l5; zeros(ng(k),1) ];
  u5 = [ u5; 
         Inf(length(ii), 1) ];
end
lc5bas = lc4Cend + 1;
lc5end = lc5bas + size(A5, 1) - 1;
% The positive reserve variables (times availability factors) are larger
% than all increment variables in all flows. 0 <= avail_fac * rPp - dPp.
% Note that these are the constraints that set
% the shadow price on the upward reserve variables
A6 = sparse(0,0); l6 = []; u6 = [];
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A6 = [ A6 ;
      sparse( [1:ng(k) 1:ng(k)]', [(rPpbas-1)+ii; (dPpbas(k):dPpend(k))' ], ...
              [avail_fac(ii, k); -ones(ng(k),1)], ng(k), nvars) ];
  l6 = [ l6 ;  zeros(ng(k),1)];
  u6 = [ u6 ; Inf(length(ii), 1) ];
end
lc6bas = lc5end + 1;
lc6end = lc6bas + size(A6,1) - 1;
% dispatches are greater than percentage of Gmax (proxy for an aggregate Pmin)
% 0 <= Pg - beta(Pc + avail_fac * rPp)
A7B = sparse(0,0); l7B = []; u7B = [];
for k = 1:nc+1
  genk = Augmented_gen(sum(ng(1:k-1)) + (1:ng(k)), :);
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  jj = find(~isload(genk) & genk(:, PMAX) > 0);
  idx = [pgbas(k):pgend(k)]';
  ngk = length(jj);
  A7B = [ A7B ;
      sparse( [  1:ngk      1:ngk           1:ngk]', ...
            [ idx(jj);      pcbas-1+ii(jj); rPpbas-1+ii(jj) ], ...
            [ ones(ngk,1);  -PminFactor(ii(jj));  -PminFactor(ii(jj)) .* avail_fac(ii(jj), k) ], ...
             ngk, nvars) ];
  l7B = [ l7B ;
         zeros(ngk,1) ];
  u7B = [ u7B;
         Inf(ngk, 1) ];
end
lc7Bbas = lc6end + 1;
lc7Bend = lc7Bbas + size(A7B,1) - 1;

% The decrement variables from the base case for generators active in kth 
% flow cannot be negative. 0 <= dPm
A8 = sparse(0,0); l8 = []; u8 = [];
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A8 = [ A8;
        sparse((1:ng(k))', (dPmbas(k):dPmend(k))', ones(ng(k),1), ng(k), nvars) ];
  l8 = [ l8; zeros(ng(k),1) ];
  u8 = [ u8; 
         Inf(length(ii), 1) ];
end
lc8bas = lc7Bend + 1;
lc8end = lc8bas + size(A8, 1) - 1;
% The negative reserve variables are larger than all decrement variables in all
% flows. 0 <= rPm - dPm . Note that these are the constraints that set
% the shadow prices on the downward reserve variables.
A9 = sparse(0,0); l9 = []; u9 = [];
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A9 = [ A9 ;
      sparse(   [ 1:ng(k)          1:ng(k)]', ...
             [(rPmbas-1)+ii; (dPmbas(k):dPmend(k))' ], ...
              [ones(ng(k),1); -ones(ng(k),1)], ng(k), nvars) ];
  l9 = [ l9 ;  zeros(ng(k),1)];
  u9 = [ u9 ; Inf(length(ii), 1) ];
end
lc9bas = lc8end + 1;
lc9end = lc9bas + size(A9,1) - 1;
% the difference between the injection and the contract
% is equal to the inc minus the dec: Pik - Pci = dPp - dPm
A10 = sparse(0,0); l10 = []; u10 = [];
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A10 = [ A10 ;
      sparse( [  1:ng(k)           1:ng(k)          1:ng(k)                 1:ng(k)]', ...
           [ (pgbas(k):pgend(k))'; pcbas(1)-1+ii;  (dPpbas(k):dPpend(k))'; (dPmbas(k):dPmend(k))'],...
           [  ones(ng(k),1);       -ones(ng(k),1); -ones(ng(k),1);          ones(ng(k),1) ], ...
              ng(k), nvars)  ];
  l10 = [ l10 ;
         zeros(ng(k),1) ];
  u10 = [ u10 ;
         zeros(ng(k),1) ];
end
lc10bas = lc9end + 1;
lc10end = lc10bas + size(A10,1) - 1;

% SECTION 5: LINEAR CONSTRAINTS FOR REACTIVE POWER RESERVES AND INCREMENTS

if HAVE_Q
  % Start defining the constraints.  First define ng(1) upward Q reserve variables
  % rqp <= rqpmax, but do it backwards so we get multipliers with
  % correct sign; MINOS writes a Lagrangian F(X) - lambda^T * G(X),
  % so we want the lower limit to be usually binding to get positive
  % multipliers, so we rewrite this as -rqpmax <= -rqp
  A11 = sparse((1:ng(1))', (rQpbas:rQpend)', -ones(ng(1),1), ng(1), nvars);
  l11 = -PositiveReactiveReserveQuantity / baseMVA;
  u11 = Inf(ng(1),1);
  lc11bas = lc10end + 1;        % linear constraint set 11 start index
  lc11end = lc11bas + ng(1)-1;  % linear constraint set 11 end index
  % Now define ng(1) downward Q reserve variables via rqm <= rqmmax,
  % or, for the reasons explained above, as -rqmmax <= -rqm.
  A12 = sparse((1:ng(1))', (rQmbas:rQmend)', -ones(ng(1),1), ng(1), nvars);
  l12 = -NegativeReactiveReserveQuantity / baseMVA;
  u12 = Inf(ng(1), 1);
  lc12bas = lc11end + 1;
  lc12end = lc12bas + ng(1) - 1;
  % The actual deviations from base flow must not exceed ramp rates
  % (for reactive, all range is allowed ....)
  A13 = sparse(0,0); l3 = []; u3 = [];
  for k = 2:nc+1
    igk = sum(ng(1:k-1)) + (1:ng(k));
    rampQk = Augmented_gen(igk, QMAX) - Augmented_gen(igk, QMIN);
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A13 = [ A13 ;
            sparse( [ 1:ng(k)           1:ng(k)]', ...
                    [ qgbas(1)-1+ii;    (qgbas(k):qgend(k))'  ], ...
                    [ -ones(ng(k),1);   ones(ng(k),1)     ], ...
                    ng(k), nvars)
          ];
    l13 = [ l13; 
            -rampQk/baseMVA ];
    u13 = [ u13;
            rampQk/baseMVA  ];
  end
  lc13bas = lc12end + 1;
  lc13end = lc13bas + size(A13,1) - 1;
  % alpha-controlled equality of Qi0 and Qc; alpha was computed earlier
  A14 = sparse([1:ng(1) 1:ng(1)]', [qgbas(1):qgend(1) qcbas:qcend]', ...
              [ones(ng(1),1); -ones(ng(1),1)],   ng(1), nvars);
  l14 = -alpha*ones(ng(1),1);
  u14 = alpha*ones(ng(1),1);
  lc14bas = lc13end + 1;
  lc14end = lc14bas + ng(1) - 1;
  % bounds on Qc
  A14B = sparse([1:ng(1)]', [qcbas:qcend]', ones(ng(1),1), ng(1), nvars);
  l14B = ReactiveContractMin / baseMVA;
  u14B = ReactiveContractMax / baseMVA;
  lc14Bbas = lc14end + 1;
  lc14Bend = lc14Bbas + ng(1) - 1;
  % The increment variables from the base case for generators active in kth 
  % flow cannot be negative. 0 <= dQp
  A15 = sparse(0,0); l15 = []; u15 = [];
  for k = 1:nc+1
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A15 = [ A15;
          sparse((1:ng(k))', (dQpbas(k):dQpend(k))', ones(ng(k),1), ng(k), nvars) ];
    l15 = [ l15; zeros(ng(k),1) ];
    u15 = [ u15; 
            Inf(length(ii), 1) ];
  end
  lc15bas = lc14Bend + 1;
  lc15end = lc15bas + size(A15, 1) - 1;
  % The positive reserve variables are larger than all increment variables in all
  % flows. 0 <= rPp - dPp. Note that these are the constraints that set
  % the shadow price on the upward reserve variables
  A16 = sparse(0,0); l16 = []; u16 = [];
  for k = 1:nc+1
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A16 = [ A16 ;
        sparse( [1:ng(k) 1:ng(k)]', [(rQpbas-1)+ii; (dQpbas(k):dQpend(k))' ], ...
                [ones(ng(k),1); -ones(ng(k),1)], ng(k), nvars) ];
    l16 = [ l16 ;  zeros(ng(k),1)];
    u16 = [ u16 ; Inf(length(ii), 1) ];
  end
  lc16bas = lc15end + 1;
  lc16end = lc16bas + size(A16,1) - 1;
  % The decrement variables from the base case for generators active in kth 
  % flow cannot be negative. 0 <= dQm
  A18 = sparse(0,0); l18 = []; u18 = [];
  for k = 1:nc+1
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A18 = [ A18;
          sparse((1:ng(k))', (dQmbas(k):dQmend(k))', ones(ng(k),1), ng(k), nvars) ];
    l18 = [ l18; zeros(ng(k),1) ];
    u18 = [ u18; 
            Inf(length(ii), 1) ];
  end
  lc18bas = lc16end + 1;
  lc18end = lc18bas + size(A18, 1) - 1;
  % The negative reserve variables are larger than all decrement variables in all
  % flows. 0 <= rQm - dQm . Note that these are the constraints that set
  % the shadow prices on the downward reserve variables.
  A19 = sparse(0,0); l19 = []; u19 = [];
  for k = 1:nc+1
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A19 = [ A19 ;
        sparse(   [ 1:ng(k)          1:ng(k)]', ...
               [(rQmbas-1)+ii; (dQmbas(k):dQmend(k))' ], ...
                [ones(ng(k),1); -ones(ng(k),1)], ng(k), nvars) ];
    l19 = [ l19 ;  zeros(ng(k),1)];
    u19 = [ u19 ; Inf(length(ii), 1) ];
  end
  lc19bas = lc18end + 1;
  lc19end = lc19bas + size(A19,1) - 1;
  % the difference between the injection and the contract
  % is equal to the inc minus the dec: Qik - Qci = dQp - dQm
  A20 = sparse(0,0); l20 = []; u20 = [];
  for k = 1:nc+1
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A20 = [ A20 ;
        sparse( [  1:ng(k)           1:ng(k)          1:ng(k)                 1:ng(k)]', ...
             [ (qgbas(k):qgend(k))'; qcbas(1)-1+ii;  (dQpbas(k):dQpend(k))'; (dQmbas(k):dQmend(k))'],...
             [ ones(ng(k),1);       -ones(ng(k),1);  -ones(ng(k),1);          ones(ng(k),1) ], ...
                ng(k), nvars)  ];
    l20 = [ l20 ;
           zeros(ng(k),1) ];
    u20 = [ u20 ;
           zeros(ng(k),1) ];
  end
  lc20bas = lc19end + 1;
  lc20end = lc20bas + size(A20,1) - 1;
else
  A11 = sparse(0,0); A12 = A11; A13 = A11; A14 = A11; A14B = A11;
  A15 = A11; A16 = A11; A18 = A11; A19 = A11; A20 = A11;
  l11 = []; l12 = []; l13 = []; l14 = []; l14B = []; l15 = [];
  l16 = []; l18 = []; l19 = []; l20 = [];
  u11 = []; u12 = []; u13 = []; u14 = []; u14B = []; u15 = [];
  u16 = []; u18 = []; u19 = []; u20 = [];
  lc11bas = lc10end + 1;   lc11end = lc11bas - 1;
  lc12bas = lc11bas;       lc12end = lc11end;
  lc13bas = lc11bas;       lc13end = lc11end;
  lc14bas = lc11bas;       lc14end = lc11end;
  lc14Bbas = lc11bas;      lc14Bend = lc11end;
  lc15bas = lc11bas;       lc15end = lc11end;
  lc16bas = lc11bas;       lc16end = lc11end;
  lc18bas = lc11bas;       lc18end = lc11end;
  lc19bas = lc11bas;       lc19end = lc11end;
  lc20bas = lc11bas;       lc20end = lc11end;
end

% total output constraints
if ~isempty(toc_map)
  A21 = sparse(0,0); l21 = []; u21 = [];
  ntoc = length(toc_cap);
  for i = 1:ntoc
    for k = 1:nc+1
      ii = find(aug_gen_stat(:, k) > 0);    % which gens active in kth flow
      if k == 1
        p = prob0;
        A21tmp = sparse(1, nvars);
      else
        jj = find(contab(:, CT_LABEL) == clist(k-1));
        p = contab(jj(1), CT_PROB);
      end
      if size(toc_map, 3) == 1      % single map for all k
        k3 = 1;
      else                          % map changes as fcn of k
        k3 = k;
      end
      kk = find(toc_map(i, ii, k3));    % which of these are in constraint i, scenario k
      jj = (pgbas(k):pgend(k))';        % active gen dispatches in kth flow
      A21tmp = A21tmp + ...
          sparse( 1, jj(kk), p * toc_map(i, ii(kk), k3)' .* toc_coeff(ii(kk), toc_type(i)), 1, nvars);
    end
    A21 = [ A21 ; A21tmp ];
  end
  l21 = -Inf(ntoc,1);
  u21 = toc_cap / baseMVA;
  lc21bas = lc20end + 1;
  lc21end = lc21bas + size(A21,1) - 1;
else
  lc21bas = lc20end + 1;    lc21end = lc20end;
  A21 = sparse(0,0); l21 = []; u21 = [];
end

% I know, stacking sparse matrices row-wise...
Acoop = [ A1; A2; A3; A4; A4B; A4C; A5; A6; A7B; A8; A9; A10; A11; A12; A13; A14; A14B; A15; A16; A18; A19; A20; A21];
lcoop = [ l1; l2; l3; l4; l4B; l4C; l5; l6; l7B; l8; l9; l10; l11; l12; l13; l14; l14B; l15; l16; l18; l19; l20; l21];
ucoop = [ u1; u2; u3; u4; u4B; u4C; u5; u6; u7B; u8; u9; u10; u11; u12; u13; u14; u14B; u15; u16; u18; u19; u20; u21];
A = [Au; Acoop];
l = [lbu; lcoop]; 
u = [ubu; ucoop];

% SECTION 6: Form generalized cost

% Now add cost on reserves; make the cost THE initial linear combination
% in MATPOWER's general cost formulation. This adds a lot of zeros to
% N, but it is easy... change to sparse N add-up of coefficients
% later for efficiency...
Cr = zeros(nvars, 1);
Cr(rPpbas:rPpend) = PositiveActiveReservePrice;
Cr(rPmbas:rPmend) = NegativeActiveReservePrice;
if HAVE_Q
  Cr(rQpbas:rQpend) = PositiveReactiveReservePrice;
  Cr(rQmbas:rQmend) = NegativeReactiveReservePrice;
end
% add costs on increments/decrements for base case
Cr(dPpbas(1):dPpend(1)) = prob0 * PositiveActiveDeltaPrice;
Cr(dPmbas(1):dPmend(1)) = prob0 * NegativeActiveDeltaPrice;
if HAVE_Q
  Cr(dQpbas(1):dQpend(1)) = prob0 * PositiveReactiveDeltaPrice;
  Cr(dQmbas(1):dQmend(1)) = prob0 * NegativeReactiveDeltaPrice;
end
if HAVE_QUADRATIC       % we have quadratic terms
  h = zeros(nvars, 1);
  h(rPpbas:rPpend) = PositiveActiveReservePrice2;
  h(rPmbas:rPmend) = NegativeActiveReservePrice2;
  if HAVE_Q
    h(rQpbas:rQpend) = PositiveReactiveReservePrice2;
    h(rQmbas:rQmend) = NegativeReactiveReservePrice2;
  end
  % add costs on increments/decrements for base case
  h(dPpbas(1):dPpend(1)) = prob0 * PositiveActiveDeltaPrice2;
  h(dPmbas(1):dPmend(1)) = prob0 * NegativeActiveDeltaPrice2;
  if HAVE_Q
    h(dQpbas(1):dQpend(1)) = prob0 * PositiveReactiveDeltaPrice2;
    h(dQmbas(1):dQmend(1)) = prob0 * NegativeReactiveDeltaPrice2;
  end
end
% now insert costs on decrements/increments for contingency flows
for k = 1:nc
  ii = find(gen_stat(:, k) > 0);
  jj = find(contab(:, CT_LABEL) == clist(k));
  p = contab(jj(1), CT_PROB);
  Cr(dPpbas(k+1):dPpend(k+1)) = p * PositiveActiveDeltaPrice(ii);
  Cr(dPmbas(k+1):dPmend(k+1)) = p * NegativeActiveDeltaPrice(ii);
  if HAVE_Q
    Cr(dQpbas(k+1):dQpend(k+1)) = p * PositiveReactiveDeltaPrice(ii);
    Cr(dQmbas(k+1):dQmend(k+1)) = p * NegativeReactiveDeltaPrice(ii);
  end
  if HAVE_QUADRATIC       % we have quadratic terms
    h(dPpbas(k+1):dPpend(k+1)) = p * PositiveActiveDeltaPrice2(ii);
    h(dPmbas(k+1):dPmend(k+1)) = p * NegativeActiveDeltaPrice2(ii);
    if HAVE_Q
      h(dQpbas(k+1):dQpend(k+1)) = p * PositiveReactiveDeltaPrice2(ii);
      h(dQmbas(k+1):dQmend(k+1)) = p * NegativeReactiveDeltaPrice2(ii);
    end
  end
end
if HAVE_QUADRATIC       % we have quadratic terms
  % let N be a sparse identity just to keep it easy to index things
  if isempty(N)
    N = speye(nvars, nvars);
    fparm = ones(nvars, 1) * [1 0 0 1];
    Cw = baseMVA * Cr;
    H = sparse(1:nvars, 1:nvars, 2 * baseMVA^2 * h, nvars, nvars);
  else
    N = [ N;  speye(nvars, nvars) ];
    fparm = [ fparm;  ones(nvars, 1) * [1 0 0 1] ];
    Cw = [ Cw;  baseMVA * Cr ];
    H = [ H sparse(size(H,1),nvars);
          sparse(nvars,size(H,2)) sparse(1:nvars, 1:nvars, 2 * baseMVA^2 * h, nvars, nvars) ];
  end
else                        % no quadratic terms
  % put the linear terms in a single row in N
  if isempty(N)
    N = baseMVA * sparse(Cr');    % make it a sparse row
    fparm = [ 1 0 0 1];           % linear identity fn, no dead zone
    Cw = 1;                       % times 1
    H  = sparse(1,1);             % no quadratic term
  else
    N = [ N;  baseMVA * sparse(Cr') ];
    fparm = [ fparm;  1 0 0 1 ];
    Cw = [ Cw;  1 ];
    H = [ H  sparse(size(H,1),1);  sparse(1,size(H,2)+1 ) ];
  end
end

% SECTION 7: Call solver
mpc = struct(...
    'baseMVA',  baseMVA, ...
    'bus',      Augmented_bus, ...
    'gen',      Augmented_gen, ...
    'branch',   Augmented_branch, ...
    'gencost',  Augmented_gencost, ...
    'A',        A, ...
    'l',        l, ...
    'u',        u, ...
    'N',        N, ...
    'fparm',    fparm, ...
    'H',        H, ...
    'Cw',       Cw ...
);
if ~isempty(iflims)
  mpc.if.map  = Augmented_ifmap;
  mpc.if.lims = Augmented_iflims;
  mpc = toggle_iflims(mpc, 'on');
end
if ~isempty(dcline)
  mpc.dcline = Augmented_dcline;
  mpc = toggle_dcline(mpc, 'on');
end
[r, success] = opf(mpc, mpopt);

[buso, geno, brancho, f, info, et] = ...
    deal(r.bus, r.gen, r.branch, r.f, r.raw.info, r.et);
if isfield(r.raw, 'dg')
    [g, jac] = deal(r.raw.g, r.raw.dg);
else
    g = [];
    jac = [];
end
if dc
  pimul = [ r.lin.mu.l.Pmis - r.lin.mu.u.Pmis;
            zeros(size(r.bus, 1), 1);
            -r.branch(:, MU_SF) * baseMVA;
            -r.branch(:, MU_ST) * baseMVA;
            r.lin.mu.l.usr  - r.lin.mu.u.usr    ];
else
  pimul = [ r.mu.nln.l      - r.mu.nln.u;
            r.lin.mu.l.usr  - r.lin.mu.u.usr ];
end

if success && (OUT_ALL > 0)
  printpf(r, 1, mpopt);
end

% SECTION 8: Preliminary unpacking of results

% create an xr that doesn't have any y variables in it
xr = [  buso(:, VA)*pi/180;
        buso(:, VM);
        geno(:, PG)/baseMVA;
        geno(:, QG)/baseMVA;
        r.var.val.z
    ];

% Pick out reserve variables from optimization vector (can be
% larger than actual reserve needed if reserve cost was zero).
% Thus result could be different from Gmax - Pc or Pc - Gmin;
% so it must be verified (see Gmin, Gmax, Qmin, Qmax below).
% For non-unity availability factors, R does not equal
% Gmax - Pc, so we have to use the optimization variables directly
PReservePlus = baseMVA * xr(rPpbas:rPpend);
PReserveMinus = baseMVA * xr(rPmbas:rPmend);
if HAVE_Q
  QReservePlus = baseMVA * xr(rQpbas:rQpend);
  QReserveMinus = baseMVA * xr(rQmbas:rQmend);
else
  QReservePlus = [];
  QReserveMinus = [];
end

% Pick out dispatches, reserve, and some multipliers, in internal order
Pik = zeros(ng(1), nc+1);
Qik = zeros(ng(1), nc+1);
lamPik = zeros(ng(1), nc+1);
lamQik = zeros(ng(1), nc+1);
dPplus = zeros(ng(1), nc+1);
dPminus = zeros(ng(1), nc+1);
Pik(:, 1) = geno(1:ng(1), PG);  % the base case dispatch
Qik(:, 1) = geno(1:ng(1), QG);
Pc = baseMVA * xr(pcbas:pcend); % the contract active dispatch
dPplus(:, 1) = baseMVA * xr(dPpbas(1):dPpend(1)); % delta P+ base case, will
dPminus(:, 1) = baseMVA * xr(dPmbas(1):dPmend(1)); % need consistency check
if HAVE_Q
  Qc = baseMVA * xr(qcbas:qcend); % the contracted reactive dispatch
  dQplus = zeros(ng(1), nc+1);
  dQminus = zeros(ng(1), nc+1);
  dQplus(:, 1) = baseMVA * xr(dQpbas(1):dQpend(1));
  dQminus(:, 1) = baseMVA * xr(dQmbas(1):dQmend(1));
end
lamPik(:, 1) = buso(geno(1:ng(1), GEN_BUS), LAM_P); % base case multipliers
lamQik(:, 1) = buso(geno(1:ng(1), GEN_BUS), LAM_Q);
lamRPp_GT_dPpik = zeros(ng(1), nc+1);
lamRPp_GT_dPpik(:, 1) = pimul(lc6bas-1+(1:ng(1))) / baseMVA;
lamRPm_GT_dPmik = zeros(ng(1), nc+1);
lamRPm_GT_dPmik(:, 1) = pimul(lc9bas-1+(1:ng(1))) / baseMVA;
temp = -pimul(lc10bas-1+(1:ng(1))) / baseMVA;
lamPikplus = zeros(ng(1), nc+1);
lamPikplus(temp > 0, 1) = temp(temp > 0);
lamPikminus = zeros(ng(1), nc+1);
lamPikminus(temp < 0, 1) = -temp(temp < 0);
lamPikp_GT_0 = zeros(ng(1), nc+1);
lamPikp_GT_0(:, 1) = pimul(lc5bas-1+(1:ng(1))) / baseMVA;
lamPikm_GT_0 = zeros(ng(1), nc+1);
lamPikm_GT_0(:, 1) = pimul(lc8bas-1+(1:ng(1))) / baseMVA;
lamPhysRamp = zeros(ng(1), nc);

%lamRPplusmax = max(0,pimul(lc1bas:lc1end)/baseMVA); % upward reserve -ramp <= -Rp+ <= 0
lamRPplusmax = pimul(lc1bas:lc1end)/baseMVA; % upward reserve -ramp <= -Rp+ <= 0
%lamRPplusmin = min(0,pimul(lc1bas:lc1end)/baseMVA); % negative means 0 limit binding
%lamRPminusmax = max(0,pimul(lc2bas:lc2end)/baseMVA); % downward reserve -ramp <= -Rp- <= 0
lamRPminusmax = pimul(lc2bas:lc2end)/baseMVA; % downward reserve -ramp <= -Rp- <= 0
%lamRPminusmin = min(0,pimul(lc2bas:lc2end)/baseMVA);
lam_alpha = -pimul(lc4bas-1+(1:ng(1))) / baseMVA;
lam_Pc = -pimul(lc4Bbas-1+(1:ng(1))) / baseMVA;

if HAVE_Q
  lamRQp_GT_dQpik = zeros(ng(1), nc);
  lamRQp_GT_dQpik(:, 1) = pimul(lc16bas-1+(1:ng(1))) / baseMVA;
  lamRQm_GT_dQmik = zeros(ng(1), nc);
  lamRQm_GT_dQmik(:, 1) = pimul(lc19bas-1+(1:ng(1))) / baseMVA;
  lamQikplus = zeros(ng(1), nc+1);
  temp = -pimul(lc20bas-1+(1:ng(1))) / baseMVA;
  lamQikplus(temp > 0, 1) = temp(temp > 0);
  lamQikminus = zeros(ng(1), nc+1);
  lamQikminus(temp < 0, 1) = -temp(temp < 0);
  lamQikp_GT_0 = zeros(ng(1), nc+1);
  lamQikp_GT_0(:, 1) = pimul(lc15bas-1+(1:ng(1))) / baseMVA;
  lamQikm_GT_0 = zeros(ng(1), nc+1);
  lamQikm_GT_0(:, 1) = pimul(lc18bas-1+(1:ng(1))) / baseMVA;
  lamPhysRampQ = zeros(ng(1), nc);
  lamRQplusmax = pimul(lc11bas:lc11end) / baseMVA;
  lamRQminusmax = pimul(lc12bas:lc12end) / baseMVA;
  lam_alphaQ = -pimul(lc14bas-1+(1:ng(1))) / baseMVA;
  lam_Qc = -pimul(lc14Bbas-1+(1:ng(1))) / baseMVA;
else
  dQplus = [];
  dQminus = [];
  lamRQp_GT_dQpik = [];
  lamRQm_GT_dQmik = [];
  lamQikplus = [];
  lamQikminus = [];
  lamQikp_GT_0 = [];
  lamQikm_GT_0 = [];
  lamPhysRampQ = [];
  lamRQplusmax = [];
  lamRQminusmax = [];
  lam_alphaQ = [];
end

% Pick out dispatches and several multipliers for each contingency.
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0);
  Pik(ii, k) = geno(sum(ng(1:k-1))+1:sum(ng(1:k)), PG); % this is how to pick from augm Gen table, works for k=1
  Qik(ii, k) = geno(sum(ng(1:k-1))+1:sum(ng(1:k)), QG); % because sum(empty) = 0.
  dPplus(ii, k) = baseMVA * xr(dPpbas(k):dPpend(k)); % delta P+ for each contingency
  dPminus(ii, k) = baseMVA * xr(dPmbas(k):dPmend(k));
  lamPik(ii, k) = buso(geno(sum(ng(1:k-1))+1:sum(ng(1:k)), GEN_BUS), LAM_P);
  lamQik(ii, k) = buso(geno(sum(ng(1:k-1))+1:sum(ng(1:k)), GEN_BUS), LAM_Q);  
  lamRPp_GT_dPpik(ii, k) = pimul((lc6bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  lamRPm_GT_dPmik(ii, k) = pimul((lc9bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  temp = -pimul((lc10bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  lamPikplus(ii( temp > 0), k) =  temp(temp > 0);
  lamPikminus(ii(temp < 0), k) = -temp(temp < 0);
  lamPikp_GT_0(ii, k) = pimul((lc5bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  lamPikm_GT_0(ii, k) = pimul((lc8bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  if k > 1
    lamPhysRamp(ii, k-1) = -pimul((lc3bas-1-ng(1)+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  end
  if HAVE_Q
    dQplus(ii, k) = baseMVA * xr(dQpbas(k):dQpend(k));
    dQminus(ii, k) = baseMVA * xr(dQmbas(k):dQmend(k));
    lamRQp_GT_dQpik(ii, k) = pimul((lc16bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
    lamRQm_GT_dQmik(ii, k) = pimul((lc19bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
    temp = -pimul((lc20bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
    lamQikplus( ii(temp > 0), k) =  temp(temp > 0);
    lamQikminus(ii(temp < 0), k) = -temp(temp < 0);
    lamQikp_GT_0(ii, k) = pimul((lc15bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
    lamQikm_GT_0(ii, k) = pimul((lc18bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
    if k > 1
      lamPhysRampQ(ii, k-1) = -pimul((lc13bas-1-ng(1)+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
    end
  end
end

% Find Gmax/Gmin, but do it only over periods where the generator
% was explicitly committed.  Thus can't do it with max()
% and min() over the matrices directly.
Gmin = zeros(ng(1), 1);
Gmax = zeros(ng(1), 1);
for i=1:ng(1)
  ii = find([ 1 gen_stat(i, :)] > 0);
  Gmin(i) = min(Pik(i, ii));
  Gmax(i) = max(Pik(i, ii));
end
Qmin = zeros(ng(1), 1);   % these are interesting even when Q is not
Qmax = zeros(ng(1), 1);   % considered in the market/cost
for i=1:ng(1)
  ii = find([ 1 gen_stat(i, :)] > 0);
  Qmin(i) = min(Qik(i, ii));
  Qmax(i) = max(Qik(i, ii));
end

% Now analyze consistency of deltas, Pc, Gmin, Gmax. Problems when
% Gmin(i) + max_k {Pikminus} + epsilon < Gmax(i) - max_k {Pikplus},
% i.e. there is a gap of epsilon between the LHS and RHS; under normal
% conditions both expressions should evaluate to Pci.
% Otherwise, if the above is tight, then Pci is equal to either
% side (with epsilon = 0).
max_dPplus = max(dPplus, [], 2);    % not sure if all are tight, but expect
max_dPminus = max(dPminus, [], 2);  % largest to be tight and equal to reserve
tight_rhs = min(Gmin + max_dPminus, ActiveContractMax);
tight_lhs = max(Gmax - max_dPplus, ActiveContractMin);
test_loose = tight_rhs-tight_lhs;
% Pc_raw = Pc;  % DEBUG CODE
for i=1:ng(1)
  if test_loose(i) > 4*baseMVA*mpopt.opf.violation
    if mpopt.verbose > 1
      fprintf('Warning: loose constraints, Pc(%d) set to ', i);
      fprintf('lower end of %.3g MW valid range\n', test_loose(i));
    end
    Pc(i) = tight_lhs(i);
  elseif test_loose(i) < -4*baseMVA*mpopt.opf.violation
    fprintf('ERROR: Constraint violation of %g MW related to inc/dec/res/Pc(%d)\n', -test_loose(i), i);
    success == 0;
%   else
%     Pc(i) = (tight_lhs(i) + tight_rhs(i))/2;
  end
  if Pc(i) > Gmax(i) && Gmax(i) > ActiveContractMax(i)
    if mpopt.verbose > 1
      fprintf('Warning: Pc(%i) > Gmax(%i) (internal order), setting to Gmax.\n', i, i);
    end
    Pc(i) = Gmax(i);
  elseif Pc(i) < Gmin(i) && Gmin(i) < ActiveContractMin(i)
    if mpopt.verbose > 1
      fprintf('Warning: Pc(%i) < Gmin(%i) (internal order), setting to Gmin.\n', i, i);
    end
    Pc(i) = Gmin(i);
  end
end
if HAVE_Q
  max_dQplus = max(dQplus, [], 2);    % not sure if all are tight, but expect
  max_dQminus = max(dQminus, [], 2);  % largest to be tight and equal to reserve
  tight_rhsq = min(Qmin + max_dQminus, ReactiveContractMax);
  tight_lhsq = max(Qmax - max_dQplus, ReactiveContractMin);
  test_looseq = tight_rhsq-tight_lhsq;
  for i=1:ng(1)
    if test_looseq(i) > 4*baseMVA*mpopt.opf.violation
      if mpopt.verbose > 1
        fprintf('Warning: loose constraints, Qc(%d) set to ', i);
        fprintf('lower end of %.3g MVAr valid range\n', test_looseq(i));
      end
      Qc(i) = tight_lhsq(i);
    elseif test_looseq(i) < -4*baseMVA*mpopt.opf.violation
      fprintf('ERROR: Constraint violation of %g MVAr related to inc/dec/res/Qc(%d)\n', -test_looseq(i), i);
      success == 0;
%     else
%       Qc(i) = (tight_lhsq(i) + tight_rhsq(i))/2;
    end
    if Qc(i) > Qmax(i) && Qmax(i) > ReactiveContractMax(i)
      if mpopt.verbose > 1
        fprintf('Warning: Qc(%i) > Qmax(%i) (in internal order), setting to Qmax.\n', i, i);
      end
      Qc(i) = Qmax(i);
    elseif Qc(i) < Qmin(i) && Qmin(i) > ReactiveContractMin(i)
      if mpopt.verbose > 1
        fprintf('Warning: Qc(%i) < Qmin(%i) (in internal order), setting to Qmin.\n', i, i);
      end
      Qc(i) = Qmin(i);
    end
  end
end
% Now Pc should be ok and can recompute dPplus, dPminus etc to make all tight
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0);
  dPplus(ii, k) = max(zeros(size(ii)), Pik(ii,k)-Pc(ii));
  dPminus(ii, k) = max(zeros(size(ii)), Pc(ii)-Pik(ii,k));
  if HAVE_Q
    dQplus(ii, k) = max(zeros(size(ii)), Qik(ii,k)-Qc(ii));
    dQminus(ii, k) = max(zeros(size(ii)), Qc(ii)-Qik(ii,k));
  end
end
% and can also recompute the reserve
% (but not with non-unity availability factors)
% PReservePlus = max(zeros(size(Pc)), Gmax - Pc);
% PReserveMinus = max(zeros(size(Pc)), Pc - Gmin);
% if HAVE_Q
%   QReservePlus = max(zeros(size(Pc)), Qmax - Qc);
%   QReserveMinus = max(zeros(size(Pc)), Qc - Qmin);
% else
%   QReservePlus = [];
%   QReserveMinus = [];
% end

% Sum some multipliers over flows
sumlamP = sum(lamPik, 2);
sumlamQ = sum(lamQik, 2);
sumlamPikplus = sum(lamPikplus, 2);
sumlamPikminus = sum(lamPikminus, 2);
sumlamRPp = sum(lamRPp_GT_dPpik, 2);
sumlamRPm = sum(lamRPm_GT_dPmik, 2);
sumlamPhysRamp = sum(lamPhysRamp, 2);
if HAVE_Q
  sumlamRQp = sum(lamRQp_GT_dQpik, 2);
  sumlamRQm = sum(lamRQm_GT_dQmik, 2);
end

% lamP0_div_prob0 = lamPik(:, 1) / prob0;


% SECTION 9: Create output results

results = struct;

% base flow results
results.base.bus = buso(1:nb(1), :);
branchwork = brancho(1:nl(1), :);
genwork = geno(inv_gen_ord, :);
genwork(:, VG) = buso(geno(inv_gen_ord, GEN_BUS), VM);
[results.base.bus, genwork, branchwork] = ...
   int2ext(i2e, results.base.bus, genwork, branchwork);
results.base.branch = branch_original;
tmp = zeros(length(base_off_branch), length([PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]));
results.base.branch(base_off_branch, [PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]) = tmp;
results.base.branch(base_on_branch, :) = branchwork;
results.base.gen = gen_original;
tmp = zeros(length(base_off_gen), length([PG;QG;MU_PMAX;MU_PMIN;MU_QMAX;MU_QMIN]));
results.base.gen(base_off_gen, [PG;QG;MU_PMAX;MU_PMIN;MU_QMAX;MU_QMIN]) = tmp;
results.base.gen(base_on_gen, :) = genwork;
if ~isempty(dcline)
  dclineo = r.dcline;
  dclinework = dclineo(1:ndc(1), :);
  dclinework(:, c.F_BUS) = i2e( dclinework(:, c.F_BUS) );
  dclinework(:, c.T_BUS) = i2e( dclinework(:, c.T_BUS) );
  results.base.dcline = dcline_original;
  tmp = zeros(length(base_off_dcline), length([c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]));
  results.base.dcline(base_off_dcline, [c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]) = tmp;
  results.base.dcline(base_on_dcline, :) = dclinework;
end
if ~isempty(iflims)
  results.base.if = iflims;
  results.base.if.P = r.if.P(1:nif(1));
  results.base.if.mu.l = r.if.mu.l(1:nif(1));
  results.base.if.mu.u = r.if.mu.u(1:nif(1));
end

% post-contingency flows results
kk = 1;  % index over actual considered contingencies
for k=1:length(clist_original) % index over original list of contingencies, 
  clabel = clist_original(k);  % including those thrown out over uncommit issues
  if any(clabel == clist) % Was this contingency actually included in analysis?
    % pull kth contingency bus data from augmented flow
    buswork = buso((1:nb(kk+1))+sum(nb(1:kk)), :);
    buswork(:, BUS_I) = (1:nb(kk+1))';
    % branch back-transformation in two steps: first to coincide with
    % "branch", which already does not have uncommitted branches on input, then
    % to coincide with branch_original, which is the original input data
    results.cont(k).branch = branch_original;
    branchwork = branch;
    branchwork1 = brancho((1:nl(kk+1))+sum(nl(1:kk)), :); %  pull from augmented
    branchwork1(:, [F_BUS;T_BUS]) = ...
        branchwork1(:, [F_BUS;T_BUS]) - sum(nb(1:kk));    % renumber FROM, TO bus #
    branchwork(branch_stat(:, kk) ~= 0, :) = branchwork1; % 1st back-transform step;
    ii = find(~branch_stat(:, kk));  % now zero out fields in inactive branches
    branchwork(ii, [BR_STATUS;PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]) = ...
      zeros(length(ii), length([BR_STATUS;PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]));
    % still have to do 2nd back-transformation but since we must first
    % reorder gens, we turn our attention to them...
    results.cont(k).gen = gen_original;
    genwork = gen;
    ii = find(gen_stat(:, kk) > 0);
    genwork(ii, :) = geno((1:ng(kk+1))+sum(ng(1:kk)), :); % pull from augm flow
    genwork(ii, GEN_BUS) = genwork(ii, GEN_BUS) - sum(nb(1:kk)); % backmap bus #
    genwork(ii, VG) = buso(geno((1:ng(kk+1))+sum(ng(1:kk))-1, GEN_BUS), VM);
    ii = find(gen_stat(:, kk) <= 0); % zero out fields for gen taken out in contingency..
    genwork(ii, [GEN_STATUS;PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]) = ...
      zeros(length(ii), length([GEN_STATUS;PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]));
    genwork = genwork(inv_gen_ord, :); % and reorder back ...
    % renumber back the buses
    [buswork, genwork, branchwork] = ...
        int2ext(i2e, buswork, genwork, branchwork);
    % and stick into results struct, zeroing out appropriate fields for
    % gens and branches that were uncommitted on input...
    results.cont(k).bus = buswork;
    results.cont(k).branch(base_on_branch, :) = branchwork;
    tmp = zeros(length(base_off_branch), length([PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]));
    results.cont(k).branch(base_off_branch, [PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]) = tmp;
    results.cont(k).gen(base_on_gen, :) = genwork;
    tmp = zeros(length(base_off_gen), length([PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]));
    results.cont(k).gen(base_off_gen, [PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]) = tmp;
    if ~isempty(dcline)
      % dcline back-transformation in two steps: first to coincide with
      % "dcline", which already does not have uncommitted dclines on input, then
      % to coincide with dcline_original, which is the original input data
      results.cont(k).dcline = dcline_original;
      dclinework = dcline;
      dclinework1 = dclineo((1:ndc(kk+1))+sum(ndc(1:kk)), :); %  pull from augmented
      dclinework1(:, [c.F_BUS;c.T_BUS]) = ...
          dclinework1(:, [c.F_BUS;c.T_BUS]) - sum(nb(1:kk));    % renumber FROM, TO bus #
%       dclinework(dcline_stat(:, kk) ~= 0, :) = dclinework1; % 1st back-transform step;
%       ii = find(~dcline_stat(:, kk));  % now zero out fields in inactive dclines
%       dclinework(ii, [c.BR_STATUS;c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]) = ...
%         zeros(length(ii), length([c.BR_STATUS;c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]));
      dclinework = dclinework1;
      % still have to do 2nd back-transformation but since we must first
      % reorder gens, we turn our attention to them...
      % renumber back the buses
      dclinework(:, c.F_BUS) = i2e( dclinework(:, c.F_BUS) );
      dclinework(:, c.T_BUS) = i2e( dclinework(:, c.T_BUS) );
      % and stick into results struct, zeroing out appropriate fields for
      % dclines that were uncommitted on input...
      results.cont(k).dcline(base_on_dcline, :) = dclinework;
      tmp = zeros(length(base_off_dcline), length([c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]));
      results.cont(k).dcline(base_off_dcline, [c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]) = tmp;
    end
    if ~isempty(iflims)
      results.cont(k).if = iflims;
      results.cont(k).if.P = r.if.P((1:nif(kk+1))+sum(nif(1:kk)));
      results.cont(k).if.mu.l = r.if.mu.l((1:nif(kk+1))+sum(nif(1:kk)));
      results.cont(k).if.mu.u = r.if.mu.u((1:nif(kk+1))+sum(nif(1:kk)));
    end
    kk = kk + 1;
  else % this contingency was deleted from the local list as a result of some
    results.cont(k).bus = [];  % equipment being uncommited on input, and
    results.cont(k).branch = []; % this contingency attempting to take it out,
    results.cont(k).gen = [];    % so we just insert empty contingency flow data
  end
end

% Compute sum of bus lambdas over contingencies plus base case
sum_bus_lam_p = buso(1:nb(1), LAM_P);
sum_bus_lam_q = buso(1:nb(1), LAM_Q);
for k = 2:nc+1
  sum_bus_lam_p = sum_bus_lam_p + buso((1:nb(k))+sum(nb(1:k-1)), LAM_P);
  sum_bus_lam_q = sum_bus_lam_q + buso((1:nb(k))+sum(nb(1:k-1)), LAM_Q);
end
results.energy.prc.sum_bus_lam_p = sum_bus_lam_p;
results.energy.prc.sum_bus_lam_q = sum_bus_lam_q;

ng_original = size(gen_original, 1);

results.energy.sum_muPmax = results.base.gen(:, MU_PMAX);
results.energy.sum_muPmin = results.base.gen(:, MU_PMIN);
for k=1:nc
  if ~isempty(results.cont(k).gen)
    results.energy.sum_muPmax = results.energy.sum_muPmax + ...
                                results.cont(k).gen(:, MU_PMAX);
    results.energy.sum_muPmin = results.energy.sum_muPmin + ...
                                results.cont(k).gen(:, MU_PMIN);
  end
end

% The following output data is all in ng_original x (nc+1) vectors; nc
% is actual number of contingencies considered in this run
tmp = zeros(ng_original, nc+1);
results.energy.delta.qty.P_pos = tmp; % incremental deltas for active power
results.energy.delta.qty.P_neg = tmp; % decremental deltas for active power
results.energy.delta.mu.P_pos_GEQ0 = tmp; % mu on inc deltas >= 0
results.energy.delta.mu.P_neg_GEQ0 = tmp;
results.reserve.mu.Rp_pos = tmp; % upward reserve for active power
results.reserve.mu.Rp_neg = tmp; % downward reserve for active power
results.energy.delta.qty.P_pos(base_on_gen, :) = dPplus(inv_gen_ord, :);
results.energy.delta.qty.P_neg(base_on_gen, :) = dPminus(inv_gen_ord, :);
results.energy.delta.mu.P_pos_GEQ0(base_on_gen, :) = lamPikp_GT_0(inv_gen_ord, :);
results.energy.delta.mu.P_neg_GEQ0(base_on_gen, :) = lamPikm_GT_0(inv_gen_ord, :);
results.reserve.mu.Rp_pos(base_on_gen, :) = lamRPp_GT_dPpik(inv_gen_ord, :);
results.reserve.mu.Rp_neg(base_on_gen, :) = lamRPm_GT_dPmik(inv_gen_ord, :);
if HAVE_Q
  results.energy.delta.qty.Q_pos = tmp;
  results.energy.delta.qty.Q_neg = tmp;
  results.energy.delta.mu.Q_pos_GEQ0 = tmp; % mu on inc deltas >= 0
  results.energy.delta.mu.Q_neg_GEQ0 = tmp;
  results.reserve.mu.Rq_pos = tmp;
  results.reserve.mu.Rq_neg = tmp;
  results.energy.delta.qty.Q_pos(base_on_gen, :) = dQplus(inv_gen_ord, :);
  results.energy.delta.qty.Q_neg(base_on_gen, :) = dQminus(inv_gen_ord, :);
  results.energy.delta.mu.Q_pos_GEQ0(base_on_gen, :) = lamQikp_GT_0(inv_gen_ord, :);
  results.energy.delta.mu.Q_neg_GEQ0(base_on_gen, :) = lamQikm_GT_0(inv_gen_ord, :);
  results.reserve.mu.Rq_pos(base_on_gen, :) = lamRQp_GT_dQpik(inv_gen_ord, :);
  results.reserve.mu.Rq_neg(base_on_gen, :) = lamRQm_GT_dQmik(inv_gen_ord, :);
end

results.energy.Pc = zeros(ng_original, 1);
results.energy.Pc(base_on_gen) = Pc(inv_gen_ord);
if HAVE_Q
  results.energy.Qc = zeros(ng_original, 1);
  results.energy.Qc(base_on_gen) = Qc(inv_gen_ord);
end
tmp = zeros(ng_original, 1);
results.energy.Gmax = tmp;
results.energy.Gmin = tmp;
results.energy.Qmax = tmp;
results.energy.Qmin = tmp;
results.energy.Gmax(base_on_gen) = Gmax(inv_gen_ord);
results.energy.Gmin(base_on_gen) = Gmin(inv_gen_ord);
results.energy.Qmax(base_on_gen) = Qmax(inv_gen_ord);
results.energy.Qmin(base_on_gen) = Qmin(inv_gen_ord);

% The following is in ng_original * nc matrices
tmp = zeros(ng_original, nc);
results.energy.mu.Ramp_P_max = tmp;
results.energy.mu.Ramp_P_max(base_on_gen, :) = lamPhysRamp(inv_gen_ord,:);
if HAVE_Q
  results.energy.mu.Ramp_Q_max = tmp;
  results.energy.mu.Ramp_Q_max(base_on_gen, :) = lamPhysRampQ(inv_gen_ord,:);
end

% The following output data is in (ng_original x 1) vectors.
tmp = zeros(ng_original, 1);
results.energy.mu.alphaP = tmp;
results.energy.mu.Pc = tmp;
results.reserve.qty.Rp_pos = tmp;
results.reserve.qty.Rp_neg = tmp;
results.reserve.prc.Rp_pos = tmp;
results.reserve.prc.Rp_neg = tmp;
results.reserve.mu.Rpmax_pos = tmp;
results.reserve.mu.Rpmax_neg = tmp;
results.energy.mu.alphaP(base_on_gen) = lam_alpha(inv_gen_ord);
results.energy.mu.Pc(base_on_gen) = lam_Pc(inv_gen_ord);
results.reserve.qty.Rp_pos(base_on_gen) = PReservePlus(inv_gen_ord);
results.reserve.qty.Rp_neg(base_on_gen) = PReserveMinus(inv_gen_ord);
results.reserve.prc.Rp_pos(base_on_gen) = sumlamRPp(inv_gen_ord);
results.reserve.prc.Rp_neg(base_on_gen) = sumlamRPm(inv_gen_ord);
results.reserve.mu.Rpmax_pos(base_on_gen) = lamRPplusmax(inv_gen_ord);
results.reserve.mu.Rpmax_neg(base_on_gen) = lamRPminusmax(inv_gen_ord);
if HAVE_Q
  results.energy.mu.alphaQ = tmp;
  results.energy.mu.Qc = tmp;
  results.reserve.qty.Rq_pos = tmp;
  results.reserve.qty.Rq_neg = tmp;
  results.reserve.prc.Rq_pos = tmp;
  results.reserve.prc.Rq_neg = tmp;
  results.reserve.mu.Rqmax_pos = tmp;
  results.reserve.mu.Rqmax_neg = tmp;
  results.energy.mu.alphaQ(base_on_gen) = lam_alphaQ(inv_gen_ord);
  results.energy.mu.Qc(base_on_gen) = lam_Qc(inv_gen_ord);
  results.reserve.qty.Rq_pos(base_on_gen) = QReservePlus(inv_gen_ord);
  results.reserve.qty.Rq_neg(base_on_gen) = QReserveMinus(inv_gen_ord);
  results.reserve.prc.Rq_pos(base_on_gen) = sumlamRQp(inv_gen_ord);
  results.reserve.prc.Rq_neg(base_on_gen) = sumlamRQm(inv_gen_ord);
  results.reserve.mu.Rqmax_pos(base_on_gen) = lamRQplusmax(inv_gen_ord);
  results.reserve.mu.Rqmax_neg(base_on_gen) = lamRQminusmax(inv_gen_ord);
end
if ~isempty(caplim_map)
    results.caplim.map = caplim_map;
    results.caplim.min = caplim_min;
    results.caplim.max = caplim_max;
    results.caplim.qty = A4C * xr * baseMVA;
    results.caplim.mu  = -pimul(lc4Cbas:lc4Cend) / baseMVA;
end
if ~isempty(toc_map)
    results.total_output.map = toc_map;
    results.total_output.cap = toc_cap;
    results.total_output.coeff = toc_coeff;
    results.total_output.type = toc_type;
    results.total_output.qty = A21 * xr * baseMVA;
    results.total_output.mu = -pimul(lc21bas:lc21end) / baseMVA;
end

% raw OPF results
results.success = success;
results.opf_results = r;

% Eval marginal cost of all injections
if HAVE_Q
  xx = [ geno(:, PG); geno(:, QG)];
else
  xx = geno(:, PG);
end
df_dPg = margcost(Augmented_gencost, xx);
if HAVE_Q
  df_dQg = df_dPg(ng_total+1:2*ng_total);
  df_dPg = df_dPg(1:ng_total);
end
%df_dPg_out = zeros(ng_original, 1);
%df_dPg_out(base_on_original) = df_dPg(inv_gen_ord);  NO

% Eval marginal cost of offers at Pc quantities
if HAVE_Q
  xx = [ Pc; Qc];
else
  xx = Pc;
end
df_dPcQc = margcost(gencost, xx);
df_dPc = zeros(ng_original, 1);
df_dPc(base_on_gen) = df_dPcQc(inv_gen_ord);
if HAVE_Q
  df_dQc = zeros(ng_original, 1);
  df_dQc(base_on_gen) = df_dPcQc(ng(1)+inv_gen_ord);
end

tmp = zeros(ng_original, 1);
sumlamP_out = tmp;
sumlamP_out(base_on_gen) = sumlamP(inv_gen_ord);
sum_genbus_lam_p_out = tmp;
sum_genbus_lam_p_out(base_on_gen) = sum_bus_lam_p(gen(inv_gen_ord, GEN_BUS));
sumlamPikplus_out = tmp;
sumlamPikplus_out(base_on_gen) = sumlamPikplus(inv_gen_ord);
sumlamPikminus_out = tmp;
sumlamPikminus_out(base_on_gen) = sumlamPikminus(inv_gen_ord);

results.energy.sumlamPikplus  = sumlamPikplus_out;
results.energy.sumlamPikminus = sumlamPikminus_out;

if success && OUT_ALL
  ngencol = 6; % number of generators per print table (max # columns)
  fprintf('\n');
  fprintf('================================================================================\n')
  fprintf('  E4ST Results\n');
  fprintf('================================================================================\n')
  fprintf('\n');
  fprintf('Total expected cost: $%g \n', f);
  fprintf('Total R+: %g \n', sum(PReservePlus));
  fprintf('Total R-: %g \n', sum(PReserveMinus));
  fprintf('Generator summary\n');
  fprintf('\n');
  for i=1:ngencol:ng_original
    igp = i:min(i+ngencol-1, ng_original);
    fprintf('           ');
    fprintf('   G%2i   ', igp);
    fprintf('\n');
    fprintf('================================================================================\n')
    fprintf('    Pc   |');
    fprintf(' %7.2f ', results.energy.Pc(igp));
    fprintf('\n');
    fprintf('   Gmin  |');
    fprintf(' %7.2f ', results.energy.Gmin(igp));
    fprintf('\n');
    fprintf('   Gmax  |');
    fprintf(' %7.2f ', results.energy.Gmax(igp));
    fprintf('\n');
%    fprintf('    P0   |');
%    fprintf(' %7.2f ', results.base.gen(igp, PG));
%    fprintf('\n');
    fprintf(' P(i, 0) |');
    fprintf(' %7.2f ', results.base.gen(igp, PG));
    fprintf('\n');
    for k = 1:nc
%      fprintf(' P(i,%2i) |', k);   % this prints k
      fprintf(' P(i,%2i) |', clist(k)); % but this prints the contingency label
      fprintf(' %7.2f ', results.cont(k).gen(igp, PG));
      fprintf('\n');
    end
    fprintf('  C''(Pc) |');
    fprintf(' %7.2f ', df_dPc(igp));
    fprintf('\n');
    fprintf(' Sum Lpg |');
    fprintf(' %7.2f ', sumlamP_out(igp));
    fprintf('\n');
    fprintf(' Lp_bus  |');
    fprintf(' %7.2f ', sum_genbus_lam_p_out(igp));
    fprintf('\n');
    fprintf('sumMuPmax|');
    fprintf(' %7.2f ', results.energy.sum_muPmax(igp));
    fprintf('\n');
    fprintf('sumMuPmin|');
    fprintf(' %7.2f ', results.energy.sum_muPmin(igp));
    fprintf('\n');
    fprintf('   Rp+   |');
    fprintf(' %7.2f ', results.reserve.qty.Rp_pos(igp));
    fprintf('\n');
    fprintf('  $Rp+   |');
    fprintf(' %7.2f ', results.reserve.prc.Rp_pos(igp));
    fprintf('\n');
    fprintf('muRp+box |');
    fprintf(' %7.2f ', results.reserve.mu.Rpmax_pos(igp));
    fprintf('\n');
    fprintf('   Rp-   |');
    fprintf(' %7.2f ', results.reserve.qty.Rp_neg(igp));
    fprintf('\n');
    fprintf('  $Rp-   |');
    fprintf(' %7.2f ', results.reserve.prc.Rp_neg(igp));
    fprintf('\n');
    fprintf('muRp-box |');
    fprintf(' %7.2f ', results.reserve.mu.Rpmax_neg(igp));
    fprintf('\n');
    fprintf('SumMuPik+|');
    fprintf(' %7.2f ', sumlamPikplus_out(igp));
    fprintf('\n');
    fprintf('SumMuPik-|');
    fprintf(' %7.2f ', sumlamPikminus_out(igp));
    fprintf('\n');
    for k = 1:nc
      if any(abs(results.energy.mu.Ramp_P_max(igp, k)) > 1e-3)
%        fprintf('muR(i,%2i)|', k);   % this prints k
        fprintf('muR(i,%2i)|', clist(k)); % but this prints the contingency label
        fprintf(' %7.2f ', results.energy.mu.Ramp_P_max(igp, k));
        fprintf('\n');
      end
    end
  end
end


%dPplus
%dPminus
%sumlamP
%lamPikplus
%sumlamPikplus
%lamPikminus
%sumlamPikminus
%lam_alpha
%lamP0_div_prob0
%sumlamQ
%lamRPp_GT_dPpik
%lamRPm_GT_dPmik
%sumlamRPp
%sumlamRPp_div_prob0 = sumlamRPp / prob0
%sumlamRPm
%sumlamRPm_div_prob0 = sumlamRPm / prob0
%lamRPplusmax
%lamRPminusmax
%lamPhysRamp
%sum_bus_lam_p
%sum_bus_lam_q


%disp('Gradient of cost wrt all injections')
%df_dPg

%disp('marginal costs at Pc/Qc quantities:')
%df_dPcQc
