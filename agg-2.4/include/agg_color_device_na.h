//----------------------------------------------------------------------------
// Anti-Grain Geometry - Version 2.4
// Copyright (C) 2002-2005 Maxim Shemanarev (http://www.antigrain.com)
// Copyright (C) 2012-2014 Sebastian Kloska (oncaphillis@snafu.de)
//
// Permission to copy, use, modify, sell and distribute this software 
// is granted provided this copyright notice appears in all copies. 
// This software is provided "as is" without express or implied
// warranty, and with no claim as to its suitability for any purpose.
//
//----------------------------------------------------------------------------
//
// Adaptation for high precision colors has been sponsored by 
// Liberty Technology Systems, Inc., visit http://lib-sys.com
//
// Liberty Technology Systems, Inc. is the provider of
// PostScript and PDF technology for software developers.
// 
//----------------------------------------------------------------------------
// Contact: mcseem@antigrain.com
//          mcseemagg@yahoo.com
//          http://www.antigrain.com
//----------------------------------------------------------------------------

#ifndef AGG_COLOR_DEVICE_NA_INCLUDED
#define AGG_COLOR_DEVICE_NA_INCLUDED

#include <math.h>
#include <stdexcept>

namespace agg
{
    // Supported byte orders for CMYK and CMYKA pixel formats
    //=======================================================================
    template<int N>
    struct order_device_n   { enum cmyk_e  { C=0, M=1, Y=2, K=3, cmyk_tag, hasAlpha=false    }; }; // ---- order_cmyk

    template<int N>
    struct order_devive_na  { enum cmyk_e  { C=0, M=1, Y=2, K=3, A=4, cmyk_tag,hasAlpha=true }; }; // ---- order_cmyka
    
    //====================================================================rgba
    template<int N>
    struct device_na
    {
        typedef double value_type;
        value_type _v[N];

        double _a; // alpha channel
        
        //--------------------------------------------------------------------
        device_na()         {
        }

        //--------------------------------------------------------------------
#if 0
        device_na(double c_, double m_, double y_, double k_,double a_=1.0) :           
            c(c_), m(m_), y(y_), k(k_),a(a_)         {
        }
#endif
   
        //--------------------------------------------------------------------
        template<int M>
        device_na(const device_na<M> & d, double a_=1.0) : _a(a_) {
            
            if(M>N)
                for(int i=0;i<N;i++)
                    _v[i]=d[i];
            else
                for(int i=0;i<M;i++)
                    _v[i]=d[i];
        }

        const value_type  operator[](int idx) const {
            if(idx>=N)
                throw std::runtime_error("illegal index");
            return _v[idx];
        }

        const value_type  alpha() const {
            return _a;
        }

        //--------------------------------------------------------------------
        void clear()        {
            for(int i=0;i<N;i++) {
                _v[i]=0.0;
            }
        }

        //--------------------------------------------------------------------

        const device_na & transparent()        {
            _a = 0.0;
            return *this;
        }

        //--------------------------------------------------------------------
        const device_na & opacity(double a_)        {
            _a = (a_ < 0.0) ? 0.0 : (a_ > 1.0)  ? 1.0  : a_;
            return *this;
        }

        //--------------------------------------------------------------------
        double opacity() const        {
            return _a;
        }

        //--------------------------------------------------------------------
        const device_na & premultiply()
        {
            for(int i=0;i<N;i++) {
                _v[i] *= _a;
            }
            return *this;
        }

        //--------------------------------------------------------------------
        const device_na &  premultiply(double a_)
        {
            if(a_ <= 0.0 || a_ <= 0.0)
            {
                clear();
                return *this;
            }

            a_ /= _a;
            for(int i=0;i<N;i++) {
                _v[i] *= a_;
            }
            _a  = a_;
            return *this;
        }

        //--------------------------------------------------------------------
        const device_na & demultiply()
        {
            if(_a == 0)
            {
                clear();
                return *this;
            }
            double a_ = 1.0 / _a;
            
            for(int i=0;i<N;i++) {
                _v[i] *= a_;
            }
            return *this;
        }
    

        //--------------------------------------------------------------------
        device_na gradient(const device_na & cl, double l) const
        {
            device_na ret;
            for(int i=0;i<N;i++) {
                ret._v[i] = _v[i] + (cl[i] - _v[i]) * l;
            }
            ret.a = _a + (cl.alpha() - _a) * l;
            return ret;
        }
#if 0        
        cmyka to_cmy() const
        {
            cmyka cl;

            cl.c = ( c * ( 1.0 - k ) + k );
            cl.m = ( m * ( 1.0 - k ) + k );
            cl.y = ( y * ( 1.0 - k ) + k );
            cl.k = 0.0;
            cl.a = a;

            return cl;
        }
        //--------------------------------------------------------------------
        static cmyka no_color() { return cmyka(0,0,0,0,0); }
        
        //--------------------------------------------------------------------
        static cmyka from_wavelength(double wl, double gamma = 1.0);
        
        //--------------------------------------------------------------------
        explicit cmyka(double wavelen, double gamma=1.0)
        {
            *this = from_wavelength(wavelen, gamma);
        }
#endif
    };

#if 0
    //----------------------------------------------------------------rgba_pre
    inline cmyka cmyk_pre(double c, double m, double y, double k,double a=1.0)
    {
        return cmyka(c, m, y, k, a).premultiply();
    }
    inline cmyka cmyk_pre(const cmyka & c)
    {
        return cmyka(c).premultiply();
    }
    inline cmyka cmyka_pre(const cmyka & c, double a)
    {
        return cmyka(c, a).premultiply();
    }
#endif
}
#endif // AGG_COLOR_DEVICE_NA_INCLUDED
