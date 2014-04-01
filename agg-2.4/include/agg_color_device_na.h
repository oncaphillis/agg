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
            
            if(M>N) {
                for(int i=0;i<N;i++)
                    _v[i]=d[i];
            }else {
                for(int i=0;i<M;i++)
                    _v[i]=d[i];

                for(int i=M;i<N;i++)
                    _v[i]=0.0;
            }
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
    }    //===================================================================cmyka8
#endif

    template<int N>
    struct device_na8
    {
        typedef int8u  value_type;
        typedef int32u calc_type;
        typedef int32  long_type;
        enum base_scale_e
        {
            base_shift = 8,
            base_scale = 1 << base_shift,
            base_mask  = base_scale - 1
        };

        typedef device_na8 self_type;

        value_type _v[N];
        value_type _a;

        //--------------------------------------------------------------------
        device_na8()
        {
        }
#if 0
        //--------------------------------------------------------------------
        cmyka8(unsigned c_, unsigned m_, unsigned y_, unsigned k_,unsigned a_=base_mask) :
            c(value_type(c_)),
            m(value_type(m_)),
            y(value_type(y_)),
            k(value_type(k_)),
            a(value_type(a_)) {}
#endif
        //--------------------------------------------------------------------
        template<int M>
        device_na8(const device_na<M> & cl, double a_) :
            _a((value_type)uround(cl.alpha()   * double(base_mask)))
        {
            if(M>N) {
                for(int i=0;i<N;i++)
                    _v[i]= uround(cl[i] * double(base_mask));
            } else {
                for(int i=0;i<M;i++)
                    _v[i]= uround(cl[i] * double(base_mask));
                for(int i=M;i<N;i++)
                    _v[i]=0;
            }
        }

        //--------------------------------------------------------------------
        template<int M>
        device_na8(const device_na8<M>& cl, unsigned a_) :
            _a(cl.alpha())
        {
            if(M>N) {
                for(int i=0;i<N;i++)
                    _v[i]=cl[i];
            } else {
                for(int i=0;i<M;i++)
                    _v[i]=cl[i];
                for(int i=M;i<N;i++)
                    _v[i]=0;
            }
        }

        //--------------------------------------------------------------------
        template<int M>
        device_na8(const device_na<M> & cl) :
            _a((value_type)uround(cl.alpha() * double(base_mask)))
        {
            if(M>N) {
                for(int i=0;i<N;i++)
                    _v[i]= uround(cl[i] * double(base_mask));
            } else {
                for(int i=0;i<M;i++)
                    _v[i]= uround(cl[i] * double(base_mask));
                for(int i=M;i<N;i++)
                    _v[i]=0;
            }
        }

        //--------------------------------------------------------------------
        void clear()
        {
            for(int i=0;i<N;i++) {
                _v[i]=0;
            }
            _a=0;
        }

        //--------------------------------------------------------------------
        const self_type& transparent()
        {
            _a = 0;
            return *this;
        }

        //--------------------------------------------------------------------
        const self_type& opacity(double a_)
        {
            _a = (a_ < 0.0) ? 0.0 : (a_ > 1.0) ? 1.0 : (value_type)uround(a_ * double(base_mask));
            return *this;
        }

        //--------------------------------------------------------------------
        double opacity() const
        {
            return double(_a) / double(base_mask);
        }

        //--------------------------------------------------------------------
        AGG_INLINE const self_type& premultiply()
        {
            if(_a == base_mask) {
                return *this;
            }

            if(_a == 0) {
                for(int j=0;j<N;j++)
                    _v[j]=0;
                return *this;
            }

            for(int j=0;j<N;j++) {
                _v[j] = value_type((calc_type(_v[j]) * _a) >> base_shift);
            }

            return *this;
        }

        //--------------------------------------------------------------------
        AGG_INLINE const self_type& premultiply(unsigned a_)
        {
            if(_a == base_mask && a_ >= base_mask)
                return *this;
            if(_a == 0 || a_ == 0)
            {
                _a = 0;
                for(int i=0;i<N;i++)
                    _v[i]=0;
                return *this;
            }

            for(int i=0;i<N;i++) {
                calc_type c_ = (calc_type(_v[i]) * a_) / _a;
                _v[i] = value_type((c_ > a_) ? a_ : c_);
            }

            _a = value_type(a_);

            return *this;
        }

        //--------------------------------------------------------------------
        AGG_INLINE const self_type& demultiply()
        {
            if(_a == base_mask) return *this;
            if(_a == 0)
            {
                for(int i=0;i<N;i++)
                    _v[i]=0;
                return *this;
            }
            for(int i=0;i<N;i++) {
                calc_type c_ = (calc_type(_v[i]) * base_mask) / _a;
                _v[i] = value_type((c_ > calc_type(base_mask)) ? calc_type(base_mask) : c_);
            }
            return *this;
        }

        //--------------------------------------------------------------------
        template<int M>
        AGG_INLINE self_type gradient(const device_na8<M>& cl, double l) const
        {
            self_type ret;
            calc_type ik = uround(l * base_scale);
            if(M>N)
            {
                for(int i=0;i<N;i++)
                {
                    ret[i] = value_type(calc_type(_v[i]) + (((calc_type(cl[i]) - _v[i]) * ik) >> base_shift));
                }
                ret.alpha(value_type(calc_type(alpha())) + (((calc_type(cl.alpha()) - alpha()) * ik) >> base_shift));
            }
            else
            {
                for(int i=0;i<M;i++)
                {
                    ret[i] = value_type(calc_type(_v[i]) + (((calc_type(cl[i]) - _v[i]) * ik) >> base_shift));
                }
                ret.alpha(value_type(calc_type(alpha())) + (((calc_type(cl.alpha()) - alpha()) * ik) >> base_shift));
                for(int i=M;i<N;i++) {
                    ret[i]=0;
                }
            }
            return ret;
        }

        //--------------------------------------------------------------------
        template<int M>
        AGG_INLINE void add(const device_na8<M>& cl, unsigned cover)
        {
            calc_type cc, ca;
            if(cover == cover_mask)
            {
                if(cl.a == base_mask)
                {
                    *this = cl;
                }
                else
                {
                    if(M>N) {
                        for(int i=0;i<N;i++)
                        {
                            cc = _v[i] + cl[i];
                            _v[i] = (cc > calc_type(base_mask)) ? calc_type(base_mask) : cc;
                        }
                        ca = _a + cl.alpha();
                        _a = (ca > calc_type(base_mask)) ? calc_type(base_mask) : ca;
                    } else {
                        for(int i=0;i<M;i++)
                        {
                            cc = _v[i] + cl[i];
                            _v[i] = (cc > calc_type(base_mask)) ? calc_type(base_mask) : cc;
                        }
                        ca = _a + cl.alpha();
                        _a = (ca > calc_type(base_mask)) ? calc_type(base_mask) : ca;
                        for(int i=M;i<N;i++)
                        {
                            _v[i] = 0;
                        }
                    }
                }
            }
            else
            {
                if(M>N) {
                    for(int i=0;i<N;i++) {
                        cc = _v[i] + ((cl[i] * cover + cover_mask / 2) >> cover_shift);
                        _v[i] = (cc > calc_type(base_mask)) ? calc_type(base_mask) : cc;
                     }
                    ca = _a + ((cl.a * cover + cover_mask / 2) >> cover_shift);
                    _a = (ca > calc_type(base_mask)) ? calc_type(base_mask) : ca;
                } else {
                    for(int i=0;i<M;i++) {
                        cc = _v[i] + ((cl[i] * cover + cover_mask / 2) >> cover_shift);
                        _v[i] = (cc > calc_type(base_mask)) ? calc_type(base_mask) : cc;
                     }

                    ca = _a + ((cl.a * cover + cover_mask / 2) >> cover_shift);
                    _a = (ca > calc_type(base_mask)) ? calc_type(base_mask) : ca;

                    for(int i=M;i<N;i++)
                    {
                        _v[i] = 0;
                    }
                }
            }
        }

        //--------------------------------------------------------------------
        template<class GammaLUT>
        AGG_INLINE void apply_gamma_dir(const GammaLUT& gamma)
        {
            for(int i=0;i<N;i++) {
                _v[i] = gamma.dir(_v[i]);
            }
        }

        //--------------------------------------------------------------------
        template<class GammaLUT>
        AGG_INLINE void apply_gamma_inv(const GammaLUT& gamma)
        {
            for(int i=0;i<N;i++) {
                _v[i] = gamma.inv(_v[i]);

            }
        }

        const value_type &alpha() const {
            return _a;
        }

        const value_type & operator[](int idx) const {
            if(idx<0 || idx>=N)
                    throw std::runtime_error("illegal idx");
            return _v[idx];
        }
        //--------------------------------------------------------------------
        static self_type no_color() { return self_type(); }
    };
#if 0
    //-------------------------------------------------------------cmyka8_pre
    inline cmyka8 cmyka8_pre(unsigned c, unsigned m, unsigned y, unsigned k,
                             unsigned a = cmyka8::base_mask)
    {
        return cmyka8(c,m,y,k,a).premultiply();
    }
    inline cmyka8 cmyka8_pre(const cmyka8 & cl)
    {
        return cmyka8(cl).premultiply();
    }
    inline cmyka8 cmyka8_pre(const cmyka8 & cl, unsigned a)
    {
        return cmyka8(cl,a).premultiply();
    }
    inline cmyka8 cmaka8_pre(const cmyka & cl)
    {
        return  cmyka8(cl).premultiply();
    }
    inline cmyka8 cmyka8_pre(const cmyka & cl, double a)
    {
        return cmyka8(cl,a).premultiply();
    }

    inline cmyka to_cmyka(const cmyka8 & cl)
    {
        return cmyka( (double) cl.c/255.0,
                      (double) cl.m/255.0,
                      (double) cl.y/255.0,
                      (double) cl.k/255.0,
                      (double) cl.a/255.0);
    }
#endif
}
#endif // AGG_COLOR_DEVICE_NA_INCLUDED
