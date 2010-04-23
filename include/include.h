//nog enkele definities:
#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#include "lapack.h"
#include "Matrix.h"
#include "BlockMatrix.h"
#include "Vector.h"
#include "BlockVector.h"
#include "TPM.h"
#include "SPM.h"
#include "PHM.h"

#include "SUP.h"
#include "EIG.h"
