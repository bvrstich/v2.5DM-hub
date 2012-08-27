//definities
//I1 Q2
#ifdef __I1Q2

#define __Q2_CON

#endif

//I1 I2
#ifdef __I

#define __I2_CON

#endif

//I1 I2 Q2
#ifdef __IQ2

#define __Q2_CON
#define __I2_CON

#endif

//I1 I2 Q1
#ifdef __IQ1

#define __I2_CON
#define __Q1_CON

#endif

//I1 I2 Q1 Q2
#ifdef __IQ

#define __I2_CON
#define __Q1_CON
#define __Q2_CON

#endif

//I1 I2 Q1 Q2 G1
#ifdef __IQG1

#define __I2_CON
#define __Q1_CON
#define __Q2_CON
#define __G1_CON

#endif

//I1 I2 Q1 Q2 G1 G2
#ifdef __IQG

#define __I2_CON
#define __Q1_CON
#define __Q2_CON
#define __G1_CON
#define __G2_CON

#endif

#include "lapack.h"
#include "Matrix.h"
#include "Vector.h"
#include "BlockMatrix.h"
#include "BlockVector.h"

#include "rxTPM.h"
#include "dDPM.h"
#include "SPM.h"
#include "TPM.h"
#include "ssdTPM.h"
#include "PHM.h"
#include "xSPM.h"
#include "dTPM.h"
#include "xTPM.h"
#include "dPPHM.h"
#include "rxPHM.h"
#include "dPHHM.h"

#include "Tools.h"

#include "SUP.h"
#include "EIG.h"
