
import pandas as pd
import numpy as np
from collections import OrderedDict
import math
from lattice_tree.LatticeTree import LatticeTree

class BondOptionPricing:
    """BondOptionPricing
    This class is used to pricing the coupon bond option include european option and American option.
    """
    """BondOptionPricing
    This class is used to pricing the coupon bond option include european option and American option.
    """
    def __init__(self,
                 T0:float,
                 T_N:float,
                 r_up_bound:float,
                 r_low_bound:float,
                 c1:float,
                 c2:float,
                 c3:float,
                 strike_price:float,
                 delta_t:float,
                 delta_T:float,
                 r:dict,
                 K:float,
                 sigma:float,
                 Rmn:dict,
                 ytm:dict
                 ):
        """
        :parameter
        :param T0: It's the maturity of the option,it's year.
        :param TN: It's the maturity of the coupon bond, it's year.
        :param r_up_bound: The r up bound in the coupon function.
        :param r_low_bound:The r low bound in the coupon function.
        :param c1:c1 in the coupon function.
        :param c2:c2 in the coupon function.
        :param c3:c3 in the coupon function.
        :param strike_price:The strike price of the option.
        :param delta_t: The delta t of the step.
        :param delta_T: The delta T of the coupon payment.
        :param r:The dictionary of the short rate r, it can be got by the LatticeTree.
        :param K:The K in the Hull-White model formula, it stands for the mean-reversion speed.
        :param sigma:The sigma in the Hull-white model formula.
        :param Rmn ytm:Rmn and ytm are the same with that in the LatticeTree and here we use Rmn and ytm to initiallize the LatticeTree.
        """
        self.T0=T0
        self.T0_step = math.floor(T0 / delta_t)
        self.T_N=T_N
        self.T_N_step=math.floor(T_N/delta_t)
        self.r_up_bound=r_up_bound
        self.r_low_bound=r_low_bound
        self.c1=c1
        self.c2=c2
        self.c3=c3
        self.strike_price=strike_price
        self.delta_t=delta_t
        self.delta_T=delta_T
        self.r=r
        self.J=math.floor(1/(2*K*delta_t))+1
        self.K=K
        self.sigma=sigma
        self.Rmn=Rmn
        self.ytm=ytm

    def f_Tj_ΔT(self,Tj_step:float):
        """
        This function is used to calculate the f(Tj;Tj,ΔT).
        :parameter
        :param Tj_step:At time Tj, what's the step.
        :return:A dictionary whose key is the state in Tj_step and value is the f(Tj;Tj,ΔT).
        """
        latticetree = LatticeTree(self.K, self.sigma, self.delta_t, self.Rmn, self.ytm)
        Tj_ΔT_step=Tj_step+math.floor(self.delta_T/self.delta_t)
        P=OrderedDict()
        P[Tj_ΔT_step] = OrderedDict()
        for i in range(-min(Tj_ΔT_step,self.J),min(Tj_ΔT_step,self.J)+1):
            P[Tj_ΔT_step][i]=1

        for j in range(Tj_ΔT_step-1,Tj_step-1,-1):
            P[j]=OrderedDict()
            for i in range(-min(j,self.J),min(j,self.J)+1):
                P[j][i] = 0
                for k in range(-min(j+1,self.J),min(j+1,self.J)+1):
                    P[j][i]+=P[j+1][k]*latticetree.get_p(i,j, k)*math.exp(-self.r[j][i]*self.delta_t)
        R=OrderedDict()
        for i in range(-min(Tj_step,self.J),min(Tj_step,self.J)+1):
            R[i]=(1/self.delta_T)*(1/(P[Tj_step][i])-1)
        return R

    def get_coupon(self)->dict:
        """
        This function is used to get the coupon.
        :return: It's a dictionary whose first key is the step and second key is the state.
        """
        r=self.r
        J=self.J
        coupon=OrderedDict()
        for j in range(self.T_N_step, -1, -1):
            coupon[j] = OrderedDict()
            R=BondOptionPricing.f_Tj_ΔT(self,j)
            if (self.T_N_step - j)%(self.delta_T/self.delta_t)>=-1e-5 and (self.T_N_step - j)%(self.delta_T/self.delta_t)<=0.5:  # Because delta_T and delta_t is not equal.
                for i in range(-min(j, J), min(j, J) + 1):
                    if R[i] >= self.r_up_bound:
                        coupon[j][i] = self.c1*self.delta_T
                    elif R[i] >= self.r_low_bound:
                        coupon[j][i] = self.c2*self.delta_T
                    elif R[i] >= 0:
                        coupon[j][i] = self.c3*self.delta_T
                    else:
                        coupon[j][i] = 0
            elif (self.T_N_step+1 - j)%(self.delta_T/self.delta_t)>=-1e-5 and (self.T_N_step+1 - j)%(self.delta_T/self.delta_t)<=0.5:
                for i in range(-min(j, J), min(j, J) + 1):
                    if R[i] >= self.r_up_bound:
                        coupon[j][i] = self.c1*self.delta_T
                    elif R[i] >= self.r_low_bound:
                        coupon[j][i] = self.c2*self.delta_T
                    elif R[i] >= 0:
                        coupon[j][i] = self.c3*self.delta_T
                    else:
                        coupon[j][i] = 0
            else:
                for i in range(-min(j, J), min(j, J) + 1):
                    coupon[j][i] = 0
        return coupon

    def get_Btc(self)->dict:
        """
        This function is used to get the Btc which is the price of the coupon bond exclude coupon, called post coupon bond.
        :return:A dictionary whose first key is step and second key is state.
        """
        J=self.J
        coupon=BondOptionPricing.get_coupon(self)
        latticetree=LatticeTree(self.K,self.sigma,self.delta_t,self.Rmn,self.ytm)
        Btc=OrderedDict()
        Btc[self.T_N_step]=OrderedDict()
        for j in range(-min(self.T_N_step,J), min(self.T_N_step,J) + 1):  # State.
            Btc[self.T_N_step][j] = 1

        for j in range(self.T_N_step-1,-1,-1):  # Step.
            Btc[j] = OrderedDict()
            for i in range(-min(j, J), min(j,J)+1):  # State.
                Btc[j][i] = 0
                for k in range(-min(j + 1, J),min(j+1,J) + 1):
                    Btc[j][i]+=(Btc[j+1][k]+coupon[j+1][k])*latticetree.get_p(i,j,k)*math.exp(-self.r[j][i]*self.delta_t)
        return Btc

    def EuropeanOptionPricing(self)->float:
        """
        This function is used to pricing the European opiton of the coupon bond.
        :return: The price of the option at time 0.
        """
        latticetree = LatticeTree(self.K, self.sigma, self.delta_t, self.Rmn, self.ytm)
        Btc=BondOptionPricing.get_Btc(self)
        V = OrderedDict()  # V is price of the option.
        V[self.T0_step] = OrderedDict()
        for i in range(-min(self.T0_step, self.J), min(self.T0_step, self.J) + 1):
            V[self.T0_step][i] = max(Btc[self.T0_step][i] - self.strike_price, 0)
        for j in range(self.T0_step - 1, -1, -1):
            V[j] = OrderedDict()
            for i in range(-min(j, self.J), min(j, self.J) + 1):
                V[j][i] = 0
                for k in range(-min(j + 1, self.J), min(j + 1, self.J) + 1):
                    V[j][i] += V[j + 1][k] * latticetree.get_p(i, j, k) * math.exp(-self.r[j][i] *self.delta_t)

        return V[0][0]

    def AmericanOPitionPricing(self)->float:
        """
        This function is used to pricing the American option of the coupon bond.
        :return: The price of the option at time 0.
        """
        latticetree = LatticeTree(self.K, self.sigma, self.delta_t, self.Rmn, self.ytm)
        Btc = BondOptionPricing.get_Btc(self)
        V = OrderedDict()  # V is price of the option.
        V[self.T0_step] = OrderedDict()
        for i in range(-min(self.T0_step, self.J), min(self.T0_step, self.J) + 1):
            V[self.T0_step][i] = max(Btc[self.T0_step][i] - self.strike_price, 0)
        for j in range(self.T0_step - 1, -1, -1):
            V[j] = OrderedDict()
            for i in range(-min(j, self.J), min(j, self.J) + 1):
                V[j][i] = 0
                for k in range(-min(j + 1, self.J), min(j + 1, self.J) + 1):
                    V[j][i] += V[j + 1][k] * latticetree.get_p(i, j, k) * math.exp(-self.r[j][i] * self.delta_t)
                V[j][i] = max(V[j][i], Btc[j][i] - self.strike_price)
        return V[0][0]


