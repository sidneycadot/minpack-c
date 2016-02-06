#include <limits>
#include <cfloat>
#include <iostream>


int main()
{
    std::cout.precision(20);

    std::cout << "*** DPMPAR values from MINPACK source ***" << std::endl;
    std::cout << std::endl;
    std::cout << "dpmpar(1) ...... : " <<  2.22044604926e-16  << std::endl;
    std::cout << "dpmpar(2) ...... : " <<  2.22507385852e-308 << std::endl;
    std::cout << "dpmpar(3) ...... : " <<  1.79769313485e308  << std::endl;
    std::cout << std::endl;

    std::cout << "*** float (C) ***" << std::endl;
    std::cout << std::endl;
    std::cout << "epsilon ...... : " << FLT_EPSILON << std::endl;
    std::cout << "min .......... : " << FLT_MIN     << std::endl;
    std::cout << "max .......... : " << FLT_MAX     << std::endl;
    std::cout << std::endl;

    std::cout << "*** double (C) ***" << std::endl;
    std::cout << std::endl;
    std::cout << "epsilon ...... : " << DBL_EPSILON << std::endl;
    std::cout << "min .......... : " << DBL_MIN     << std::endl;
    std::cout << "max .......... : " << DBL_MAX     << std::endl;
    std::cout << std::endl;

    std::cout << "*** long double (C) ***" << std::endl;
    std::cout << std::endl;
    std::cout << "epsilon ...... : " << LDBL_EPSILON << std::endl;
    std::cout << "min .......... : " << LDBL_MIN     << std::endl;
    std::cout << "max .......... : " << LDBL_MAX     << std::endl;
    std::cout << std::endl;

    std::cout << "*** float (C++) ***" << std::endl;
    std::cout << std::endl;
    std::cout << "epsilon ...... : " << std::numeric_limits<float>::epsilon() << std::endl;
    std::cout << "min .......... : " << std::numeric_limits<float>::min() << std::endl;
    std::cout << "max .......... : " << std::numeric_limits<float>::max() << std::endl;
    std::cout << std::endl;

    std::cout << "*** double (C++) ***" << std::endl;
    std::cout << std::endl;
    std::cout << "epsilon ...... : " << std::numeric_limits<double>::epsilon() << std::endl;
    std::cout << "min .......... : " << std::numeric_limits<double>::min() << std::endl;
    std::cout << "max .......... : " << std::numeric_limits<double>::max() << std::endl;
    std::cout << std::endl;

    std::cout << "*** long double (C++) ***" << std::endl;
    std::cout << std::endl;
    std::cout << "epsilon ...... : " << std::numeric_limits<long double>::epsilon() << std::endl;
    std::cout << "min .......... : " << std::numeric_limits<long double>::min() << std::endl;
    std::cout << "max .......... : " << std::numeric_limits<long double>::max() << std::endl;
    std::cout << std::endl;

    return 0;
}
