
import pandas as pd
import numpy as np
from collections import OrderedDict
import math
from bootstrapping.Bootstrapping import Bootstrapping
from lattice_tree.LatticeTree import LatticeTree
from option_pricing.BondOptionPricing import BondOptionPricing
from montecarlo.MonteCarlo import MonteCarlo
import warnings
warnings.filterwarnings("ignore")

def main():
    Term = [2, 3, 4, 5, 7, 10, 20, 30]
    swap_rate = [0.00202, 0.00288, 0.00434, 0.00597, 0.00919, 0.01263, 0.017, 0.01788]
    Rmn = dict(zip(Term, swap_rate))
    ytm = OrderedDict()
    ytm[1] = 0.00206
    delta_t=1
    delta_T=0.5
    # --------------------------------------------------------------
    # boostrap
    # --------------------------------------------------------------
    # bootstrap=Bootstrapping(Rmn,ytm)
    # ytm_discount_data=bootstrap.get_ytm_discount_data()
    # bootstrap.draw_yield_curve()
    # bootstrap.draw_discount_curve()

    # --------------------------------------------------------------
    # lattiec tree
    # --------------------------------------------------------------
    latticetree=LatticeTree(0.12,0.035,delta_t,Rmn,ytm)
    # r=latticetree.get_r_theta_dict()[0]
    # theta=latticetree.get_theta_data()
    # theta.to_excel(r'./theta.xlsx')
    latticetree.draw_r_picture()
    latticetree.draw_theta_picture()

    # --------------------------------------------------------------
    # bond option pricing by lattice tree
    # --------------------------------------------------------------
    # bond_option_pricing = BondOptionPricing(T0=0.5, T_N=2, r_up_bound=100, r_low_bound=-100, c1=0.1, c2=0.1,
    #                                         c3=0.1, strike_price=0.2, delta_t=delta_t, delta_T=delta_T, r=r,
    #                                         K=0.12, sigma=0.035, Rmn=Rmn, ytm=ytm)
    # ep=bond_option_pricing.EuropeanOptionPricing()
    # ap=bond_option_pricing.AmericanOPitionPricing()
    # print('The price of the European option is ',ep)
    # print('The price of the American option is ', ap)

    # --------------------------------------------------------------
    # bond option pricing by monte carlo
    # --------------------------------------------------------------
    # montecarlo=MonteCarlo(ytm_discount_data,5,15,0.5)
    # df_montecarlo=pd.DataFrame()
    # price=montecarlo.run_monte_carlo(2000)
    # print(price)

if __name__=="__main__":
    main()