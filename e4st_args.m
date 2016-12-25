function [baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
        Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q, opf_results, ...
        energy, avail_fac, toc_map, toc_cap, toc_coeff, toc_type, dr2aap_idx, dr2aap_factor] = ...
    e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, ...
        Au, lbu, ubu, mpopt, N, fparm, H, Cw)
%E4ST_ARGS  Parses and initializes input arguments for E4ST_SOLVE.
%   [BASEMVA, BUS, GEN, BRANCH, GENCOST, OFFER, CONTAB, AU, LBU, UBU, ...
%       MPOPT, N, FPARM, H, CW, HAVE_Q, OPF_RESULTS, ENERGY] = E4ST_ARGS(...)
%   Returns the full set of initialized input arguments for E4ST_SOLVE,
%   filling in default values for missing arguments. See Examples below for
%   the possible calling syntax options.
%
%   Examples:
%       Output argument options:
%
%       [baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
%           Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q] = e4st_args(...)
%       [baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
%           Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q, ...
%           opf_results, energy] = e4st_args(...)
%       [baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
%           Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q, ...
%           opf_results, energy, avail_fac] = e4st_args(...)
%       [baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
%           Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q, ...
%           maxbuild_map, maxbuild_capacity] = e4st_args(...)
%       [baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
%           Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q, ...
%           maxbuild_map, maxbuild_capacity, avail_fac] = e4st_args(...)
%       [baseMVA, bus, gen, branch, gencost, dcline, iflims, offer, contab, ...
%           Au, lbu, ubu, mpopt, N, fparm, H, Cw, HAVE_Q, ...
%           maxbuild_map, maxbuild_capacity, avail_fac, ...
%           toc_map, toc_cap, toc_coeff, toc_type] = e4st_args(...)
%
%       Input argument options:
%
%       e4st_args(mpc)
%       e4st_args(mpc, mpopt)
%       e4st_args(mpc, offer, contab)
%       e4st_args(mpc, offer, contab, mpopt)
%       e4st_args(mpc, offer, contab, A, l, u)
%       e4st_args(mpc, offer, contab, A, l, u, mpopt)
%       e4st_args(mpc, offer, contab, A, l, u, mpopt, N, fparm, H, Cw)
%       e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab)
%       e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, ...
%                                   mpopt)
%       e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, ...
%                                   A, l, u)
%       e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, ...
%                                   A, l, u, mpopt)
%       e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, ...
%                                   A, l, u, mpopt, N, fparm, H, Cw)
%
%   mpc is a MATPOWER case file or case struct with the fields baseMVA, bus,
%   gen, branch, gencost, and (optionally) areas. It may also include a
%   'contingencies' field (in place of contab argument). The offer argument
%   can be a struct or a matrix. If it is a struct, it has the following fields
%   for active power quantities, each an ng x 1 vector ...
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
%   Values in the offer and contab arguments override any corresponding
%   values in mpc.
%
%   The outputs are the MATPOWER case data, an offer matrix, a contab,
%   user-defined constraints, MATPOWER options struct, user-defined
%   costs, a flag indicating whether the offers include any reactive
%   power portion, and the 'opf_results' and 'energy' fields from the
%   results of a previous E4ST_SOLVE run. If the mpc struct passed in
%   includes a 'maxbuild' field (with sub-fields 'map' and 'capacity'),
%   then they are returned in the fields named 'maxbuild_map' and
%   'maxbuild_capacity'. If the mpc struct passed in includes an
%   'availability_factor' field, (an ng x 1 vector or ng x (nc+1) matrix)
%   then it is returned as 'avail_fac'. If the input mpc struct includes
%   an 'total_output' field, with sub-fields 'map', 'cap', 'coeff' and
%   (optionally) 'type' then they are returned in the arguments 'toc_map',
%   'toc_cap', 'toc_coeff' and 'toc_type', respectively.

%   E4ST
%   Copyright (c) 2000-2016 by Power System Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://e4st.com/ for more info.


if ischar(baseMVA) || isstruct(baseMVA)  %% passing filename or struct
  %---- e4st_args(baseMVA, bus,   gen,    branch, areas, gencost, offer, contab, Au,    lbu, ubu, mpopt, N, fparm, H, Cw)
  % 11  e4st_args(mpc,     offer, contab, A,      l,     u,       mpopt, N,      fparm, H,   Cw)
  %  7  e4st_args(mpc,     offer, contab, A,      l,     u,       mpopt)
  %  6  e4st_args(mpc,     offer, contab, A,      l,     u)
  %  4  e4st_args(mpc,     offer, contab, mpopt)
  %  3  e4st_args(mpc,     offer, contab)
  %  2  e4st_args(mpc,     mpopt)
  %  1  e4st_args(mpc)
  if any(nargin == [1, 2, 3, 4, 6, 7])
    casefile = baseMVA;
    if nargin == 11
      Cw    = ubu;
      H     = lbu;
      fparm = Au;
      N     = contab;
      mpopt = offer;
      ubu   = gencost;
      lbu   = areas;
      Au    = branch;
      contab = gen;
      offer = bus;
    elseif nargin == 7
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      mpopt = offer;
      ubu   = gencost;
      lbu   = areas;
      Au    = branch;
      contab = gen;
      offer = bus;
    elseif nargin == 5
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      mpopt = mpoption;
      ubu   = gencost;
      lbu   = areas;
      Au    = branch;
      contab = gen;
      offer = bus;
    elseif nargin == 4
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      ubu   = []; lbu = []; Au = sparse(0, 0);
      mpopt = branch;
      contab = gen;
      offer = bus;
    elseif nargin == 3
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      ubu   = []; lbu = []; Au = sparse(0, 0);
      mpopt = mpoption;
      contab = gen;
      offer = bus;
    elseif nargin == 2
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      ubu   = []; lbu = []; Au = sparse(0, 0);
      contab = [];
      offer = [];
      mpopt = bus;
    elseif nargin == 1
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      ubu   = []; lbu = []; Au = sparse(0, 0);
      contab = [];
      offer = [];
      mpopt = mpoption;
    end
  else
    error('e4st_args.m: Incorrect input parameter order, number or type');
  end
  mpc = loadcase(casefile);
  [baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
  if isfield(mpc, 'availability_factor')
    avail_fac = mpc.availability_factor;
  else
    avail_fac = [];
  end
  if isfield(mpc, 'dcline') && toggle_dcline(mpc, 'status')
    dcline = mpc.dcline;
  else
    dcline = [];
  end
  if isfield(mpc, 'if') && toggle_iflims(mpc, 'status')
    iflims = mpc.if;
  else
    iflims = [];
  end
  if ~isempty(avail_fac)
    if size(avail_fac, 1) ~= size(gen, 1)
      error('e4st_args.m: availability_factor must have the same number of rows as gen');
    end
    if any(any(avail_fac < 0 | avail_fac > 1))
      error('e4st_args.m: all elements x of availability_factor must satisfy 0 <= x <= 1');
    end
  end
  if isfield(mpc, 'maxbuild') && isfield(mpc.maxbuild, 'capacity') && ...
        isfield(mpc.maxbuild, 'map')
    opf_results = mpc.maxbuild.map;
    energy      = mpc.maxbuild.capacity;
    if size(opf_results, 2) ~= size(gen, 1)
      error('e4st_args.m: number of maxbuild.map cols must equal number of gens');
    end
    if size(opf_results, 1) ~= size(energy, 1)
      error('e4st_args.m: maxbuild.map and maxbuild.capacity must have the same number of rows');
    end
  else
    if isfield(mpc, 'opf_results')
      opf_results = mpc.opf_results;
    else
      opf_results = [];
    end
    if isfield(mpc, 'energy')
      energy = mpc.energy;
    else
      energy = [];
    end
  end
  if isfield(mpc, 'total_output') && ...
        isfield(mpc.total_output, 'map') && ...
        isfield(mpc.total_output, 'cap') && ...
        isfield(mpc.total_output, 'coeff') && ...
        ~isempty(mpc.total_output.map) && ...
        ~isempty(mpc.total_output.cap) && ...
        ~isempty(mpc.total_output.coeff);
    toc_map     = mpc.total_output.map;
    toc_cap     = mpc.total_output.cap;
    toc_coeff   = mpc.total_output.coeff;
    if isfield(mpc.total_output, 'type')
      toc_type = mpc.total_output.type;
    else
      toc_type = ones(size(toc_map, 1), 1);
    end
    if size(toc_map, 2) ~= size(gen, 1)
      error('e4st_args.m: number of total_output.map cols must equal number of gens');
    end
    if size(toc_map, 1) ~= size(toc_cap, 1)
      error('e4st_args.m: total_output.map and total_output.cap must have the same number of rows');
    end
    if size(toc_coeff, 1) ~= size(gen, 1)
      error('e4st_args.m: total_output.coeff must have same number of rows as gen');
    end
    if any(toc_type) < 0 || any(toc_type) > size(toc_coeff, 2) || any(toc_type ~= fix(toc_type))
      error('e4st_args.m: total_output.type be integer indices corresponding to columns of total_output.coeff');
    end
  else
    toc_map     = [];
    toc_cap     = [];
    toc_coeff   = [];
    toc_type    = [];
  end

  % Apply demand response to annual average price
  if isfield(mpc, 'dr2aap') && ...
        isfield(mpc.dr2aap, 'load_idx') && ...
        ~isempty(mpc.dr2aap.load_idx) && ...
        isfield(mpc.dr2aap, 'load_factor') && ...
        ~isempty(mpc.dr2aap.load_factor)
      dr2aap_idx = mpc.dr2aap.load_idx;
      dr2aap_factor = mpc.dr2aap.load_factor;
  else
    dr2aap_idx = [];
    dr2aap_factor = [];
  end

  if isempty(contab)
    if isfield(mpc, 'contingencies')   %% take contingencies from mpc
      contab = mpc.contingencies;
    else
      error('e4st_args.m: no ''contab'' specified');
    end
  end
  if isempty(offer)
    offer = e4st_offer2mat(mpc);
  end
else    %% passing individual data matrices
  %---- e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, Au,   lbu, ubu, mpopt, N, fparm, H, Cw)
  % 16  e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, Au,   lbu, ubu, mpopt, N, fparm, H, Cw)
  % 12  e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, Au,   lbu, ubu, mpopt)
  % 11  e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, Au,   lbu, ubu)
  %  9  e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab, mpopt)
  %  8  e4st_args(baseMVA, bus, gen, branch, areas, gencost, offer, contab)
  dcline    = [];
  iflims    = [];
  opf_results = [];
  energy    = [];
  avail_fac = [];
  toc_map   = [];
  toc_cap   = [];
  toc_coeff = [];
  toc_type  = [];
  if any(nargin == [8, 9, 11, 12])
    if nargin == 12
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
    elseif nargin == 11
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      mpopt = mpoption;
    elseif nargin == 9
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      mpopt = Au;
      ubu   = []; lbu = []; Au = sparse(0, 0);
    elseif nargin == 8
      Cw    = []; H = sparse(0, 0); fparm = []; N = sparse(0, 0);
      ubu   = []; lbu = []; Au = sparse(0, 0);
      mpopt = mpoption;
    end
  else
    error('e4st_args.m: Incorrect input parameter order, number or type');
  end
end

%% put offer into matrix until we are done with gen re-ordering
ng = size(gen, 1);
offer = e4st_offer2mat(offer, ng);

%% do we have any reactive power offers?
if size(offer, 1) == 2*ng
  rr = (ng+1):(2*ng);
  if any(any( offer(rr, 1:10) )) || ...
        any( offer(rr, 11) ~= -Inf ) || ...
        any( offer(rr, 12) ~= Inf )
    HAVE_Q = 1;
  else              %% all zero offers, unbounded Qc
    offer = offer(1:ng, :);
    HAVE_Q = 0;
  end
else
  HAVE_Q = 0;
end
