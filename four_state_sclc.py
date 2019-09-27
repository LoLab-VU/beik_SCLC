# @authors sbeik, lharris

from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

Monomer('NE')
Monomer('NEv1')
Monomer('NEv2')
Monomer('NonNE')

Parameter('NE_init', 100)
Initial(NE(), NE_init)

Observable('NE_obs', NE())
Observable('NEv1_obs', NEv1())
Observable('NEv2_obs', NEv2())
Observable('NonNE_obs', NonNE())

# NE and effects
Parameter('k_ne_div', 1)
Parameter('k_ne_die', 0.9)
Parameter('ne_plus_nonNE_div_percent', 0.01)
Parameter('ne_plus_nonNE_die_percent', 0.01)
Expression('k_ne_plus_nonNe_div', ne_plus_nonNE_div_percent*(k_ne_div)*NonNE_obs)
Expression('k_ne_plus_nonNe_die', (ne_plus_nonNE_die_percent*(k_ne_die)*NonNE_obs))
Expression('k_ne_div_tot',(k_ne_div+k_ne_plus_nonNe_div))
Expression('k_ne_die_tot',(k_ne_die-k_ne_plus_nonNe_die))
Rule('NE_div', NE() >> NE() + NE(), k_ne_div_tot)
Rule('NE_die', NE() >> None, k_ne_die_tot)

# NEv1 and effects
Parameter('k_nev1_div', 1)
Parameter('k_nev1_die', 0.9)
Parameter('nev1_plus_nonNE_div_percent',0.01)
Parameter('nev1_plus_nonNE_die_percent',0.01)
Expression('k_nev1_plus_nonNe_div', nev1_plus_nonNE_div_percent*(k_nev1_div)*NonNE_obs)
Expression('k_nev1_plus_nonNe_die', (nev1_plus_nonNE_die_percent*(k_nev1_die)*NonNE_obs))
Expression('k_nev1_div_tot',(k_nev1_div+k_nev1_plus_nonNe_div))
Expression('k_nev1_die_tot',(k_nev1_die+k_nev1_plus_nonNe_die))
Rule('NEv1_div', NEv1() >> NEv1() + NEv1(), k_nev1_div_tot)
Rule('NEv1_die', NEv1() >> None, k_nev1_die_tot)

# NEv2 and effects
Parameter('k_nev2_div', 1)
Parameter('k_nev2_die', 0.9)
Parameter('nev2_plus_nonNE_div_percent',0.01)
Parameter('nev2_plus_nonNE_die_percent',0.01)
Expression('k_nev2_plus_nonNe_div', nev2_plus_nonNE_div_percent*(k_nev2_div)*NonNE_obs)
Expression('k_nev2_plus_nonNe_die', (nev2_plus_nonNE_die_percent*(k_nev2_die)*NonNE_obs))
Expression('k_nev2_div_tot',k_nev2_div+k_nev2_plus_nonNe_div)
Expression('k_nev2_die_tot',(k_nev2_die-k_nev2_plus_nonNe_die))
Rule('NEv2_div', NEv2() >> NEv2() + NEv2(), k_nev2_div_tot)
Rule('NEv2_die', NEv2() >> None, k_nev2_die_tot)

# nonNE and effects
Parameter('k_nonNe_div', 0.9)
Parameter('k_nonNe_die', 1.0)
Parameter('ne_decr_nonNE_div_percent',0.000)
Parameter('nev1_decr_nonNE_div_percent',0.000)
Parameter('nev2_decr_nonNE_div_percent',0.000)
Expression('k_ne_decr_nonNe_div', ne_decr_nonNE_div_percent*(k_nonNe_div)*NE_obs)
Expression('k_nev1_decr_nonNe_div', nev1_decr_nonNE_div_percent*(k_nonNe_div)*NEv1_obs)
Expression('k_nev2_decr_nonNe_div', nev2_decr_nonNE_div_percent*(k_nonNe_div)*NEv2_obs)
Expression('k_nonNe_div_tot',(k_nonNe_div-(k_ne_decr_nonNe_div+k_nev1_decr_nonNe_div+k_nev2_decr_nonNe_div)))
Rule('NonNE_div', NonNE() >> NonNE() + NonNE(), k_nonNe_div_tot)
Rule('NonNE_die', NonNE() >> None, k_nonNe_die)

# NE to NEv1 transition and effects
Parameter('kf_diff_ne_nev1', 0.1)
Parameter('kr_diff_ne_nev1', 0.1)
Parameter('ne_plus_nonNE_diff_nev1_percent', 0.01)
Expression('kf_diff_ne_nev1_plus_nonNe', ne_plus_nonNE_diff_nev1_percent*(kf_diff_ne_nev1)*NonNE_obs)
Expression('kf_diff_ne_nev1_tot',(kf_diff_ne_nev1+kf_diff_ne_nev1_plus_nonNe))
Rule('NE_diff_NEv1', NE() | NEv1(), kf_diff_ne_nev1_tot, kr_diff_ne_nev1)

# NE to NEv2 transition and effects
Parameter('kf_diff_ne_nev2', 0.1)
Parameter('kr_diff_ne_nev2', 0.1)
Parameter('ne_plus_nonNE_diff_nev2_percent',0.01)
Expression('kf_diff_ne_nev2_plus_nonNe', ne_plus_nonNE_diff_nev2_percent*(kf_diff_ne_nev2)*NonNE_obs)
Expression('kf_diff_ne_nev2_tot',(kf_diff_ne_nev2+kf_diff_ne_nev2_plus_nonNe))
Rule('NE_diff_NEv2', NE() | NEv2(), kf_diff_ne_nev2_tot, kr_diff_ne_nev2)

# NEv1 to NEv2 transition (no effects - no current evidence that another subtype helps or hinders this)
Parameter('kf_diff_nev1_nev2', 0.1)
Parameter('kr_diff_nev1_nev2', 0.1)
Rule('NEv1_diff_NEv2', NEv1() | NEv2(), kf_diff_nev1_nev2, kr_diff_nev1_nev2)

# NEv1 to nonNE transition (no effects - no current evidence that another subtype helps or hinders this)
Parameter('kf_diff_nev1_nonNe', 0.1)
Rule('NEv1_diff_NonNE', NEv1() >> NonNE(), kf_diff_nev1_nonNe)

tspan = np.linspace(0, 100, 101)

sim = ScipyOdeSimulator(model, verbose=False)
x = sim.run(tspan)

plt.figure()

for obs in model.observables:
    plt.plot(tspan, x.all[obs.name], label=obs.name)

plt.xlabel('time')
plt.ylabel('cell count')
plt.legend(loc=0)
plt.tight_layout()

cell_tot = np.array([sum(x.observables[i]) for i in range(len(x.observables))])

plt.figure()
plt.fill_between(tspan, x.all[model.observables[0].name] / cell_tot, label=model.observables[0].name)

sum_prev = x.all[model.observables[0].name]
for i in range(1,len(model.observables)-1):
    plt.fill_between(tspan, (x.all[model.observables[i].name] + sum_prev) / cell_tot, sum_prev / cell_tot, label=model.observables[i].name)
    sum_prev += x.all[model.observables[i].name]


plt.fill_between(tspan, [1]*len(tspan), sum_prev / cell_tot, label=model.observables[-1].name)
plt.xlabel('time')
plt.ylabel('cell fraction')
plt.legend(loc=0)
plt.tight_layout()

plt.show()
