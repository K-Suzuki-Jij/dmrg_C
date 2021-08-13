//
//  FIND_SITE_SZ.c
//  1DXXZ_VF_DMRG
//
//  Created by Kohei Suzuki on 2019/06/13.
//  Copyright © 2019 Kohei Suzuki. All rights reserved.
//

#include "Header.h"

int FIND_SITE_SZ(int basis, int spin) {
   
   return spin - 2*basis;
  
}
