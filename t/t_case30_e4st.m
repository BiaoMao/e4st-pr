function mpc = t_case30_e4st
%T_CASE30_E4ST   Modified case30.m for security constrained OPF.
%   Please see CASEFORMAT for details on the case file format.
%
% This is based on case30V19.m with the following modifications:
%   1. Generator costs (not including dispatchable loads) have
%      been converted back to single segment piece-wise linear
%      in order to test the automatic convert to polynomial
%      done by OPF.
%   ... plus additional changes.

%   E4ST

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	135	1	1.05	0.95;
	2	2	0	0	0	0	1	1	0	135	1	1.1	0.95;
	3	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	4	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	5	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	6	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	7	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	8	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	9	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	10	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	11	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	12	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	13	2	0	0	0	0	2	1	0	135	1	1.1	0.95;
	14	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	15	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	16	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	17	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	18	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	19	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	20	1	0	0	0	0	2	1	0	135	1	1.05	0.95;
	21	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	22	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	23	2	0	0	0	0	2	1	0	135	1	1.1	0.95;
	24	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	25	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	26	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	27	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	28	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	29	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	30	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	5.88	0	15	-10	1	100	1	20	0	0	20	-14	18	-6	8	20	20	20	170	0;
	1	10.3	0	26.25	-17.5	1	100	1	35	0	0	35	-24.5	31.5	-10.5	14	35	35	35	170	0;
	2	19.05	0	18.75	-12.5	1	100	1	25	0	0	25	-17.5	22.5	-7.5	10	25	25	25	80	0;
	2	22.87	0	22.5	-15	1	100	1	30	0	0	30	-21	27	-9	12	30	30	30	80	0;
	22	12.95	0	22.5	-15	1	100	1	30	0	0	30	-21	27	-9	12	3	3	3	77.5	0;
	22	8.64	0	15	-10	1	100	1	20	0	0	20	-14	18	-6	8	20	20	20	77.5	0;
	27	17.12	0	26.25	-17.5	1	100	1	35	0	0	35	-24.5	31.5	-10.5	14	3.5	3.5	3.5	63.7	0;
	27	9.79	0	15	-10	1	100	1	20	0	0	20	-14	18	-6	8	20	20	20	63.7	0;
	23	10.47	0	15	-10	1	100	1	20	0	0	20	-14	18	-6	8	2	2	2	50	0;
	23	18.33	0	26.25	-17.5	1	100	1	35	0	0	35	-24.5	31.5	-10.5	14	7	7	7	50	0;
	13	25.62	0	22.5	-15	1	100	1	30	0	0	30	-21	27	-9	12	3	3	3	59.7	0;
	13	29.88	0	26.25	-17.5	1	100	1	35	0	0	35	-24.5	31.5	-10.5	14	7	7	7	59.7	0;
% dispatchable loads
	2	-21.7	-12.7	0	-12.7	1	100	1	0	-21.7	0	0	0	0	0	0	100	100	100	0	0;
	3	-2.4	-1.2	0	-1.2	1	100	1	0	-2.4	0	0	0	0	0	0	100	100	100	0	0;
	4	-7.6	-1.6	0	-1.6	1	100	1	0	-7.6	0	0	0	0	0	0	100	100	100	0	0;
	7	-22.8	-10.9	0	-10.9	1	100	1	0	-22.8	0	0	0	0	0	0	100	100	100	0	0;
	8	-30	-15	0	-15	1	100	1	0	-30	0	0	0	0	0	0	100	100	100	0	0;
	10	-5.8	-2	0	-2	1	100	1	0	-5.8	0	0	0	0	0	0	100	100	100	0	0;
	12	-11.2	-7.5	0	-7.5	1	100	1	0	-11.2	0	0	0	0	0	0	100	100	100	0	0;
	14	-6.2	-1.6	0	-1.6	1	100	1	0	-6.2	0	0	0	0	0	0	100	100	100	0	0;
	15	-8.2	-2.5	0	-2.5	1	100	1	0	-8.2	0	0	0	0	0	0	100	100	100	0	0;
	16	-3.5	-1.8	0	-1.8	1	100	1	0	-3.5	0	0	0	0	0	0	100	100	100	0	0;
	17	-9	-4.8	0	-4.8	1	100	1	0	-9	0	0	0	0	0	0	100	100	100	0	0;
	18	-3.2	-0.9	0	-0.9	1	100	1	0	-3.2	0	0	0	0	0	0	100	100	100	0	0;
	19	-9.5	-3.4	0	-3.4	1	100	1	0	-9.5	0	0	0	0	0	0	100	100	100	0	0;
	20	-2.2	-0.7	0	-0.7	1	100	1	0	-2.2	0	0	0	0	0	0	100	100	100	0	0;
	21	-17.5	-11.2	0	-11.2	1	100	1	0	-17.5	0	0	0	0	0	0	100	100	100	0	0;
	23	-3.2	-1.6	0	-1.6	1	100	1	0	-3.2	0	0	0	0	0	0	100	100	100	0	0;
	24	-8.7	-6.7	0	-6.7	1	100	1	0	-8.7	0	0	0	0	0	0	100	100	100	0	0;
	26	-3.5	-2.3	0	-2.3	1	100	1	0	-3.5	0	0	0	0	0	0	100	100	100	0	0;
	29	-2.4	-0.9	0	-0.9	1	100	1	0	-2.4	0	0	0	0	0	0	100	100	100	0	0;
	30	-10.6	-1.9	0	-1.9	1	100	1	0	-10.6	0	0	0	0	0	0	100	100	100	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.02	0.06	0.03	130	130	130	1	0	1	-360	360;
	1	3	0.05	0.19	0.02	130	130	130	1	0	1	-360	360;
	2	4	0.06	0.17	0.02	65	65	65	1	0	1	-360	360;
	3	4	0.01	0.04	0	130	130	130	1	0	1	-360	360;
	2	5	0.05	0.2	0.02	130	130	130	1	0	1	-360	360;
	2	6	0.06	0.18	0.02	65	65	65	1	0	1	-360	360;
	4	6	0.01	0.04	0	90	90	90	1	0	1	-360	360;
	5	7	0.05	0.12	0.01	70	70	70	1	0	1	-360	360;
	6	7	0.03	0.08	0.01	130	130	130	1	0	1	-360	360;
	6	8	0.005	0.02	0	66	66	66	1	0	1	-360	360;
	6	9	0	0.21	0	65	65	65	1	0	1	-360	360;
	6	10	0	0.56	0	20	20	20	1	0	1	-360	360;
	9	11	0	0.21	0	65	65	65	1	0	1	-360	360;
	9	10	0	0.11	0	30	30	30	1	0	1	-360	360;
	4	12	0	0.13	0	45	45	45	1	0	1	-360	360;
	12	13	0	0.14	0	90	90	90	1	0	1	-360	360;
	12	14	0.12	0.26	0	32	32	32	1	0	1	-360	360;
	12	15	0.07	0.13	0	32	32	32	1	0	1	-360	360;
	12	16	0.09	0.2	0	32	32	32	1	0	1	-360	360;
	14	15	0.22	0.2	0	16	16	16	1	0	1	-360	360;
	16	17	0.08	0.19	0	16	16	16	1	0	1	-360	360;
	15	18	0.11	0.22	0	16	16	16	1	0	1	-360	360;
	18	19	0.06	0.13	0	16	16	16	1	0	1	-360	360;
	19	20	0.03	0.07	0	32	32	32	1	0	1	-360	360;
	10	20	0.09	0.21	0	17	17	17	1	0	1	-360	360;
	10	17	0.03	0.08	0	32	32	32	1	0	0	-360	360;
	10	21	0.03	0.07	0	32	32	32	1	0	1	-360	360;
	10	22	0.07	0.15	0	32	32	32	1	0	1	-360	360;
	21	22	0.005	0.01	0	64	64	64	1	0	1	-360	360;
	15	23	0.05	0.1	0	45	45	45	1	0	1	-360	360;
	22	24	0.12	0.18	0	16	16	16	1	0	1	-360	360;
	23	24	0.13	0.27	0	20	20	20	1	0	1	-360	360;
	24	25	0.19	0.33	0	16	16	16	1	0	1	-360	360;
	25	26	0.25	0.38	0	16	16	16	1	0	1	-360	360;
	25	27	0.11	0.21	0	16	16	16	1	0	1	-360	360;
	28	27	0	0.4	0	25	25	25	1	0	1	-360	360;
	27	29	0.22	0.42	0	16	16	16	1	0	1	-360	360;
	27	30	0.32	0.6	0	16	16	16	1	0	1	-360	360;
	29	30	0.24	0.45	0	16	16	16	1	0	1	-360	360;
	8	28	0.06	0.2	0.02	32	32	32	1	0	1	-360	360;
	6	28	0.02	0.06	0.01	32	32	32	1	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	1	0	0	2	0	0	1	80;
	1	0	0	2	0	0	1	95;
	1	0	0	2	0	0	1	80;
	1	0	0	2	0	0	1	95;
	1	0	0	2	0	0	1	5;
	1	0	0	2	0	0	1	55;
	1	0	0	2	0	0	1	5;
	1	0	0	2	0	0	1	55;
	1	0	0	2	0	0	1	5;
	1	0	0	2	0	0	1	25;
	1	0	0	2	0	0	1	5;
	1	0	0	2	0	0	1	25;
% dispatchable loads
	2	0	0	2	10000	0	0	0;
	2	0	0	2	10000	0	0	0;
	2	0	0	2	10000	0	0	0;
	2	0	0	2	10000	0	0	0;
	2	0	0	2	10000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
	2	0	0	2	5000	0	0	0;
];

%%-----  DC Line Data  -----
%	fbus	tbus	status	Pf	Pt	Qf	Qt	Vf	Vt	Pmin	Pmax	QminF	QmaxF	QminT	QmaxT	loss0	loss1
mpc.dcline = [
	1	23	1	10	8.9	0	0	1	1	1	10	-10	10	-10	10	0.9	0.01;
	2	13	1	0	0	0	0	1	1	0	1	0	0	0	0	0	0.01;
	22	27	1	0	0	0	0	1	1	1	10	-10	10	-10	10	0	0.01;
];

%%-----  Interface Flow Limit Data  -----%%
%	ifnum	branchidx (negative defines opposite direction)
mpc.if.map = [
	1	-12;	%% 1 : area 1 imports
	1	-14;
	1	-15;
	1	-36;
	2	15;	%% 2 : area 2 imports
	2	25;
	2	26;
	2	-32;
	3	12;	%% 3 : area 3 imports
	3	14;
	3	-25;
	3	-26;
	3	32;
	3	36;
];

%% DC model flow limits in MW
%% (negative and positive directions can be different)
%	ifnum	lower	upper
% mpc.if.lims = [
% 	1	-15	84.4;	%% area 1 imports
% 	2	-59.7	20;	%% area 2 imports
% ];
mpc.if.lims = [
	1	-15	85;	%% area 1 imports
	2	-60	20;	%% area 2 imports
];

%%-----  Emission Constraint Data  -----%%
ng = size(mpc.gen, 1);
mpc.total_output.map = zeros(3, ng);
mpc.total_output.map(1, 1:4)  = 1;      %% area 1
mpc.total_output.map(2, 9:12) = 1;      %% area 2
mpc.total_output.map(3, 5:8)  = 1;      %% area 3
mpc.total_output.cap = [100; 90; 9.5];
mpc.total_output.coeff = zeros(ng, 2);
mpc.total_output.coeff(1:12, 1) = 1;
mpc.total_output.coeff(1:12, 2) = 0.1;
mpc.total_output.type = [1; 1; 2];
