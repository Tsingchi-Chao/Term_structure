
import pandas as pd
import numpy as np
from collections import OrderedDict
import math
from bootstrapping.Bootstrapping import Bootstrapping
from lattice_tree.LatticeTree import LatticeTree
from option_pricing.BondOptionPricing import BondOptionPricing

def main():
    Term = [2, 3, 4, 5, 7, 10, 20, 30]
    swap_rate = [0.00202, 0.00288, 0.00434, 0.00597, 0.00919, 0.01263, 0.017, 0.01788]
    Rmn = dict(zip(Term, swap_rate))
    ytm = OrderedDict()
    ytm[1] = 0.00206
    #--------------------------------------------------------------
    # analysis of delta t on option price
    #--------------------------------------------------------------
    # data_European_optoin=pd.DataFrame()
    # data_American_option=pd.DataFrame()
    # delta_T=1
    # for delta_t in np.arange(0.5,0.05,-0.01):
    #     print(delta_t)
    #     latticetree = LatticeTree(0.12, 0.035, delta_t, Rmn, ytm)
    #     r = latticetree.get_r_theta_dict()[0]
    #     bond_option_pricing=BondOptionPricing(T0=0.5,T_N=2,r_up_bound=0.2,r_low_bound=0.05,c1=0,c2=0,
    #                                           c3=0,strike_price=0.2,delta_t=delta_t,delta_T=delta_T,r=r,
    #                                           K=0.12,sigma=0.035,Rmn=Rmn,ytm=ytm)
    #     ep = bond_option_pricing.EuropeanOptionPricing()
    #     ap = bond_option_pricing.AmericanOPitionPricing()
    #     data_European_optoin.loc[delta_t,'European_option_price_adj']=ep
    #     data_American_option.loc[delta_t,'American_option_price_adj']=ap
    #     print(str(delta_t)+"对应的欧式期权价格:"+str(ep))
    #     print(str(delta_t)+"对应的美式期权价格“"+str(ap))
    # data_European_optoin.to_excel(r'./European_optoin_delta_t.xlsx')
    # data_American_option.to_excel(r'./American_option_delta_t.xlsx')

    # --------------------------------------------------------------
    # analysis of delta T on option price
    # --------------------------------------------------------------
    # data_European_optoin = pd.DataFrame()
    # data_American_option = pd.DataFrame()
    # delta_t = 0.05
    # latticetree = LatticeTree(0.12, 0.035, delta_t, Rmn, ytm)
    # r = latticetree.get_r_theta_dict()[0]
    # for multiplier in np.arange(1,40):
    #     print(multiplier)
    #     delta_T = multiplier * delta_t
    #     bond_option_pricing=BondOptionPricing(T0=0.5,T_N=2,r_up_bound=0.2,r_low_bound=0.05,c1=0.1,c2=0.2,
    #                                           c3=0.3,strike_price=0.2,delta_t=delta_t,delta_T=delta_T,r=r,
    #                                           K=0.12,sigma=0.035,Rmn=Rmn,ytm=ytm)
    #     ep = bond_option_pricing.EuropeanOptionPricing()
    #     ap = bond_option_pricing.AmericanOPitionPricing()
    #     data_European_optoin.loc[delta_T, 'European_option_price'] = ep
    #     data_American_option.loc[delta_T, 'American_option_price'] = ap
    # data_European_optoin.to_excel(r'./European_optoin_delta_big_T.xlsx')
    # data_American_option.to_excel(r'./American_option_delta_big_T.xlsx')

    # --------------------------------------------------------------
    # analysis of strike price on option price
    # --------------------------------------------------------------
    # data_European_optoin = pd.DataFrame()
    # data_American_option = pd.DataFrame()
    # delta_t = 0.1
    # delta_T=0.2
    # latticetree = LatticeTree(0.12, 0.035, delta_t, Rmn, ytm)
    # r = latticetree.get_r_theta_dict()[0]
    # for strike_price in np.arange(0.05,0.5,0.01):
    #     print(strike_price)
    #     bond_option_pricing=BondOptionPricing(T0=0.5,T_N=2,r_up_bound=0.2,r_low_bound=0.05,c1=0.1,c2=0.2,
    #                                           c3=0.3,strike_price=strike_price,delta_t=delta_t,delta_T=delta_T,r=r,
    #                                           K=0.12,sigma=0.035,Rmn=Rmn,ytm=ytm)
    #     ep = bond_option_pricing.EuropeanOptionPricing()
    #     ap = bond_option_pricing.AmericanOPitionPricing()
    #     data_European_optoin.loc[strike_price, 'European_option_price'] = ep
    #     data_American_option.loc[strike_price, 'American_option_price'] = ap
    # data_European_optoin.to_excel(r'./European_optoin_delta_strike_price.xlsx')
    # data_American_option.to_excel(r'./American_option_delta_strike_price.xlsx')

    # --------------------------------------------------------------
    # analysis of r_up_bound price on option price
    # --------------------------------------------------------------
    # data_European_optoin = pd.DataFrame()
    # data_American_option = pd.DataFrame()
    # delta_t = 0.1
    # delta_T=0.2
    # latticetree = LatticeTree(0.12, 0.035, delta_t, Rmn, ytm)
    # r = latticetree.get_r_theta_dict()[0]
    # for r_up_bound in np.arange(0.06,0.5,0.01):
    #     print(r_up_bound)
    #     bond_option_pricing=BondOptionPricing(T0=0.5,T_N=2,r_up_bound=r_up_bound,r_low_bound=0.05,c1=0.1,c2=0.2,
    #                                           c3=0.3,strike_price=0.2,delta_t=delta_t,delta_T=delta_T,r=r,
    #                                           K=0.12,sigma=0.035,Rmn=Rmn,ytm=ytm)
    #     ep = bond_option_pricing.EuropeanOptionPricing()
    #     ap = bond_option_pricing.AmericanOPitionPricing()
    #     data_European_optoin.loc[r_up_bound, 'European_option_price'] = ep
    #     data_American_option.loc[r_up_bound, 'American_option_price'] = ap
    # data_European_optoin.to_excel(r'./European_optoin_delta_r_up_bound.xlsx')
    # data_American_option.to_excel(r'./American_option_delta_r_up_bound.xlsx')

    # --------------------------------------------------------------
    # analysis of r_low_bound price on option price
    # --------------------------------------------------------------
    # data_European_optoin = pd.DataFrame()
    # data_American_option = pd.DataFrame()
    # delta_t = 0.1
    # delta_T=0.2
    # latticetree = LatticeTree(0.12, 0.035, delta_t, Rmn, ytm)
    # r = latticetree.get_r_theta_dict()[0]
    # for r_low_bound in np.arange(0.01,0.19,0.005):
    #     print(r_low_bound)
    #     bond_option_pricing=BondOptionPricing(T0=0.5,T_N=2,r_up_bound=0.2,r_low_bound=r_low_bound,c1=0.1,c2=0.2,
    #                                           c3=0.3,strike_price=0.2,delta_t=delta_t,delta_T=delta_T,r=r,
    #                                           K=0.12,sigma=0.035,Rmn=Rmn,ytm=ytm)
    #     ep = bond_option_pricing.EuropeanOptionPricing()
    #     ap = bond_option_pricing.AmericanOPitionPricing()
    #     data_European_optoin.loc[r_low_bound, 'European_option_price'] = ep
    #     data_American_option.loc[r_low_bound, 'American_option_price'] = ap
    # data_European_optoin.to_excel(r'./European_optoin_delta_r_low_bound.xlsx')
    # data_American_option.to_excel(r'./American_option_delta_r_low_bound.xlsx')

    # --------------------------------------------------------------
    # analysis of T_N price on option price
    # --------------------------------------------------------------
    # data_European_optoin = pd.DataFrame()
    # data_American_option = pd.DataFrame()
    # delta_t = 0.1
    # delta_T=0.2
    # latticetree = LatticeTree(0.12, 0.035, delta_t, Rmn, ytm)
    # r = latticetree.get_r_theta_dict()[0]
    # for T_N in np.arange(0.6,3.5,0.1):
    #     print(T_N)
    #     bond_option_pricing=BondOptionPricing(T0=0.5,T_N=T_N,r_up_bound=0.2,r_low_bound=0.05,c1=0.1,c2=0.2,
    #                                           c3=0.3,strike_price=0.2,delta_t=delta_t,delta_T=delta_T,r=r,
    #                                           K=0.12,sigma=0.035,Rmn=Rmn,ytm=ytm)
    #     ep = bond_option_pricing.EuropeanOptionPricing()
    #     ap = bond_option_pricing.AmericanOPitionPricing()
    #     data_European_optoin.loc[T_N, 'European_option_price'] = ep
    #     data_American_option.loc[T_N, 'American_option_price'] = ap
    # data_European_optoin.to_excel(r'./European_option_delta_T_N.xlsx')
    # data_American_option.to_excel(r'./American_option_delta_T_N.xlsx')

    # --------------------------------------------------------------
    # analysis of T0 price on option price
    # --------------------------------------------------------------
    data_European_optoin = pd.DataFrame()
    data_American_option = pd.DataFrame()
    delta_t = 0.1
    delta_T=0.2
    latticetree = LatticeTree(0.12, 0.035, delta_t, Rmn, ytm)
    r = latticetree.get_r_theta_dict()[0]
    for T0 in np.arange(0.1,2,0.1):
        print(T0)
        bond_option_pricing=BondOptionPricing(T0=T0,T_N=2,r_up_bound=0.2,r_low_bound=0.05,c1=0.1,c2=0.2,
                                              c3=0.3,strike_price=0.2,delta_t=delta_t,delta_T=delta_T,r=r,
                                              K=0.12,sigma=0.035,Rmn=Rmn,ytm=ytm)
        ep = bond_option_pricing.EuropeanOptionPricing()
        ap = bond_option_pricing.AmericanOPitionPricing()
        data_European_optoin.loc[T0, 'European_option_price'] = ep
        data_American_option.loc[T0, 'American_option_price'] = ap
    data_European_optoin.to_excel(r'./European_option_delta_T_0.xlsx')
    data_American_option.to_excel(r'./American_option_delta_T_0.xlsx')


if __name__=="__main__":
    main()





