# Term_structure
这部分是term structure的课程作业,有一些重要的代码可以后面继续使用和参考。

## Pricing_option_of_Bond
该部分在Ho-Lee模型下，采用蒙特卡洛模拟的方法，为Bond进行了定价。

## Project
该project是利率期限结构的大作业，其中分为四个部分。

* bootstrapping
  即根据市场信息根据线性插值法获得zero-coupon bond yield curve
  
* lattice_tree
  根据Hull-white short rate model，得到其lattice tree。
  
 * option_pricing
   根据得到的lattice tree为期权进行定价
   
 * montecarlo
   同时利用蒙特卡洛进行定价，与lattice tree定价结果可以进行对比。
