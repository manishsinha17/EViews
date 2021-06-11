smpl @all 'the whole sample
genr y=100*dlog(adj_close) 'generate return series, 1/3/2014--1/30/2015
scalar alpha=0.10 ' set VaR level at 5%
scalar q=@qnorm(alpha) 'inverse CDF of standard normal at alpha

smpl 1/4/2012 1/31/2018 'first rolling sample
scalar T=@obssmpl 'rolling sample size

smpl 02/01/2018 4/20/2021 'out-of-sample period
scalar N=@obssmpl 'out-of-sample size

genr VaR_GARCH=0 'to store VaR from GARCH etc.
genr VaR_TGARCH=0
genr VaR_EGARCH=0

genr Qhat_GARCH=0 'to store the empirical check function at each out-of-sample point from GARCH etc.
genr Qhat_TGARCH=0
genr Qhat_EGARCH=0

genr alphahat_GARCH=0 'to store the empirical coverage at each out-of-sample point from GARCH etc.
genr alphahat_TGARCH=0
genr alphahat_EGARCH=0

genr hfor_GARCH=0 'to store one-step-ahead volatility from GARCH etc.
genr hfor_TGARCH=0
genr hfor_EGARCH=0

for !i=1 to N  'use "for" loop
smpl @first+!i @first+!i+T-1 'i-th rolling sample, note that @first is N/A for y
genr yy=y 'the in-sample data is passed to yy
'the following chunk is to estimate GARCH and forecast VaR
smpl @first+!i @first+!i+T-1 'rolling sample
equation eq1.arch(1,1) yy c 'estimate ARMA(0,)+GARCH(1,1) with the equation named "eq1"
eq1.makegarch h 'create fitted volatility from GARCH
genr uhat=resid 'create residuals
smpl @first+!i+T @first+!i+T 'one-day-ahead out-of-sample period
hfor_GARCH=@coefs(2)+@coefs(3)*uhat(-1)^2+@coefs(4)*h(-1) 'generate volatility forecast from GARCH
VaR_GARCH=@coefs(1)+(hfor_GARCH^0.5)*q 'generate one-day-ahead VaR

'the following chunk is to estimate TGARCH and forecast VaR
smpl @first+!i @first+!i+T-1
equation eq2.arch(1,1,t) yy c  'the option "t" is for TGARCH; use "e" for EGARCH
eq2.makegarch h 'create fitted volatility from TGARCH
genr uhat=resid
smpl @first+!i+T @first+!i+T
hfor_TGARCH=@coefs(2)+@coefs(3)*uhat(-1)^2+@coefs(4)*uhat(-1)^2*(uhat(-1)<0)+@coefs(5)*h(-1) 'generate volatility forecast from TGARCH
VaR_TGARCH=@coefs(1)+(hfor_TGARCH^0.5)*q

'the following chunk is to estimate EGARCH and forecast VaR
smpl @first+!i @first+!i+T-1
equation eq3.arch(1,1,e) yy c  'the option "t" is for TGARCH; use "e" for EGARCH
eq3.makegarch h 'create fitted volatility from EGARCH
genr uhat=resid
smpl @first+!i+T @first+!i+T
hfor_EGARCH=exp(@coefs(2)+@coefs(3)*abs(uhat(-1)/@SQRT(h(-1)))+@coefs(4)*uhat(-1)/@SQRT(h(-1))+@coefs(5)*log(h(-1))) 'generate volatility forecast from TGARCH
VaR_EGARCH=@coefs(1)+(hfor_EGARCH^0.5)*q
next

smpl 02/01/2018 4/20/2021  'out-of-sample period
group VaR VaR_GARCH VaR_TGARCH VaR_EGARCH  'group VaR forecasts from GARCH and TAGARCH and name the group "VaR"
group Qhat Qhat_GARCH Qhat_TGARCH Qhat_EGARCH 'group the check function Qhat (now 0) from VaR_GARCH and VaR_TAGARCH and name the group "Qhat"
group alphahat alphahat_GARCH alphahat_TGARCH  alphahat_EGARCH 'group alphahat (now 0) from VaR_GARCH and VaR_TAGARCH and name the group "alphahat"
for !i=1 to 3 'assign values to Qhat and  alphahat
alphahat(!i)=y<VaR(!i)
Qhat(!i)=(alpha-(y<VaR(!i)))*(y-VaR(!i))
next

freeze(QhatStatTable) Qhat.stats 'show descriptive statistics; the mean gives the empirical out-of-sample check function value; freeze as a table
freeze(alphahatStatTable) alphahat.stats 'show descriptive statistics; the mean gives the empirical out-of-sample coverage probability; freeze as a table

group y_VaR y VaR
freeze(yVaR) y_VaR.line 'plot VaR and y together for the out-of-sample period; freeze as a figure

delete yy uhat h eq1 eq2 eq3 open close low high volume 'delete objects not needed

Then do these : 

genr garchE=uhat2-hfor_garch
genr tgarchE=uhat2-hfor_tgarch
genr garche2=garche^2
genr tgarche2=tgarche^2
genr egarche=uhat2-hfor_egarch
genr egarche2=egarche^2
















