# Introduction #

we update the isotopic peaks difference parameters by using the experiment data of ISB\_control dataset.

In the commonddef.h, modified the code:
//define the isotopic mass difference
const static double ISO\_DIFF[MAX\_ISO](MAX_ISO.md)={
> /**0.0,statics parameters**/
  1. 002844793399012,
> 2.005150790031021,
> 3.005776092238567,
> 4.004639208298311,
> 4.999369511325370,
> //original
> /**1.003308559119778,
> 2.0062822183398778,
> 3.007692024454673,
> 4.009119148582318,
> 5.010832505631107**/
> //therory
> //0.0,
> //1.0033548,//mass of neutron
> //2.0067096,
> //3.0100644,
> //4.0134192,
> //5.0167740
};