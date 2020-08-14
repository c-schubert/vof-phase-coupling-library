/*
This file contains the source term functions for lorentz forces and joule heat

License (MIT):

Copyright (c) 2016-2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/ 
#include "vof_pc_main.h"


  DEFINE_SOURCE(Jouleheating, c, t, dS, eqn)
  {
    #if (N_UDM >= 0)
    return C_UDMI(c, t, UDM_JH);
    #else
    return 0
    #endif
  }


DEFINE_SOURCE(v_x_lorentz, c, t, dS, eqn)
{
  #if (N_UDM >= 1)
  return C_UDMI(c, t, UDM_LFx);
  #else
  return 0;
  #endif
}


DEFINE_SOURCE(v_y_lorentz, c, t, dS, eqn)
{
  #if (N_UDM >= 2)
  return C_UDMI(c, t, UDM_LFy);
  #else
  return 0;
  #endif
}


DEFINE_SOURCE(v_z_lorentz, c, t, dS, eqn)
{
  #if (N_UDM >= 3)
  return C_UDMI(c, t, UDM_LFz);
  #else
  return 0;
  #endif
}

