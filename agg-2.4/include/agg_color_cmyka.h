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

#ifndef AGG_COLOR_CMYKA_INCLUDED
#define AGG_COLOR_CMYKA_INCLUDED

#include <math.h>

namespace agg
{
    // Supported byte orders for CMYK and CMYKA pixel formats
    //=======================================================================
    
    struct order_cmyk   { enum cmyk_e  { C=0, M=1, Y=2, K=3, cmyk_tag, hasAlpha=false    }; }; // ---- order_cmyk
    struct order_cmyka  { enum cmyk_e  { C=0, M=1, Y=2, K=3, A=4, cmyk_tag,hasAlpha=true }; }; // ---- order_cmyka
    
    //====================================================================rgba

    struct cmyka
    {
        typedef double value_type;

        double c;
        double m;
        double y;
        double k;
        double a;
        
        //--------------------------------------------------------------------
        cmyka() 
        {
        }

        //--------------------------------------------------------------------
        cmyka(double c_, double m_, double y_, double k_,double a_=1.0) :           
            c(c_), m(m_), y(y_), k(k_),a(a_) 
        {
        }
        
        //--------------------------------------------------------------------
        cmyka(const cmyka & c, double a_) : 
            c(c.c), m(c.m), y(c.y), k(c.k), a(a_) 
        {
        }

        //--------------------------------------------------------------------
        void clear()
        {
            c = m = y = k = a = 0;
        }

        //--------------------------------------------------------------------
        const cmyka & transparent()
        {
            a = 0.0;
            return *this;
        }

        //--------------------------------------------------------------------
        const cmyka & opacity(double a_)
        {
            if(a_ < 0.0) 
                a_ = 0.0;
            if(a_ > 1.0) 
                a_ = 1.0;
            a = a_;
            return *this;
        }

        //--------------------------------------------------------------------
        double opacity() const
        {
            return a;
        }

        //--------------------------------------------------------------------
        const cmyka & premultiply()
        {
            c *= a;
            m *= a;
            y *= a;
            k *= a;
            return *this;
        }

        //--------------------------------------------------------------------
        const cmyka &  premultiply(double a_)
        {
            if(a <= 0.0 || a_ <= 0.0)
            {
                c = m = y = k = a = 0.0;
                return *this;
            }
            a_ /= a;
            c *= a_;
            m *= a_;
            y *= a_;
            k *= a_;

            a  = a_;
            return *this;
        }

        //--------------------------------------------------------------------
        const cmyka & demultiply()
        {
            if(a == 0)
            {
                c = m = y = k = 0;
                return *this;
            }
            double a_ = 1.0 / a;
            c *= a_;
            m *= a_;
            y *= a_;
            k *= a_;
            return *this;
        }
    

        //--------------------------------------------------------------------
        cmyka gradient(cmyka cl, double l) const
        {
            cmyka ret;
            ret.c = c + (cl.c - c) * l;
            ret.m = m + (cl.m - m) * l;
            ret.y = y + (cl.y - y) * l;
            ret.k = k + (cl.k - k) * l;
            ret.a = a + (cl.a - a) * l;
            return ret;
        }
        
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
    };

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

    //===================================================================cmyka8

    struct cmyka8
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
        typedef cmyka8 self_type;


        value_type c;
        value_type m;
        value_type y;
        value_type k;
        value_type a;

        //--------------------------------------------------------------------
        cmyka8() 
        {
        }

        //--------------------------------------------------------------------
        cmyka8(unsigned c_, unsigned m_, unsigned y_, unsigned k_,unsigned a_=base_mask) :
            c(value_type(c_)), 
            m(value_type(m_)), 
            y(value_type(y_)), 
            k(value_type(k_)), 
            a(value_type(a_)) {}

        //--------------------------------------------------------------------

        cmyka8(const cmyka & cl, double a_) :
            c((value_type)uround(cl.c * double(base_mask))), 
            m((value_type)uround(cl.m * double(base_mask))), 
            y((value_type)uround(cl.y * double(base_mask))), 
            k((value_type)uround(cl.k * double(base_mask))), 
            a((value_type)uround(a_   * double(base_mask))) 
        {
        }

        //--------------------------------------------------------------------
        cmyka8(const self_type& cl, unsigned a_) :
            c(cl.c), 
            m(cl.m), 
            y(cl.y), 
            k(cl.k),
            a(value_type(a_)) 
        {
        }

        //--------------------------------------------------------------------
        cmyka8(const cmyka & cl) :
            c((value_type)uround(cl.c * double(base_mask))), 
            m((value_type)uround(cl.m * double(base_mask))), 
            y((value_type)uround(cl.y * double(base_mask))), 
            k((value_type)uround(cl.k * double(base_mask))), 
            a((value_type)uround(cl.a * double(base_mask))) {}

        //--------------------------------------------------------------------
        void clear()
        {
            c = m = y = k = a = 0;
        }
        
        //--------------------------------------------------------------------
        const self_type& transparent()
        {
            a = 0;
            return *this;
        }

        //--------------------------------------------------------------------
        const self_type& opacity(double a_)
        {
            if(a_ < 0.0) a_ = 0.0;
            if(a_ > 1.0) a_ = 1.0;
            a = (value_type)uround(a_ * double(base_mask));
            return *this;
        }

        //--------------------------------------------------------------------
        double opacity() const
        {
            return double(a) / double(base_mask);
        }

        //--------------------------------------------------------------------
        AGG_INLINE const self_type& premultiply()
        {
            if(a == base_mask) 
            {
                return *this;
            }
            if(a == 0)
            {
                c = m = y = k = 0;
                return *this;
            }

            c = value_type((calc_type(c) * a) >> base_shift);
            m = value_type((calc_type(m) * a) >> base_shift);
            y = value_type((calc_type(y) * a) >> base_shift);
            k = value_type((calc_type(k) * a) >> base_shift);

            return *this;
        }

        //--------------------------------------------------------------------
        AGG_INLINE const self_type& premultiply(unsigned a_)
        {
            if(a == base_mask && a_ >= base_mask) 
                return *this;
            if(a == 0 || a_ == 0)
            {
                c = m = y = k = a = 0;
                return *this;
            }

            calc_type c_ = (calc_type(c) * a_) / a;
            calc_type m_ = (calc_type(m) * a_) / a;
            calc_type y_ = (calc_type(y) * a_) / a;
            calc_type k_ = (calc_type(k) * a_) / a;

            c = value_type((c_ > a_) ? a_ : c_);
            m = value_type((m_ > a_) ? a_ : m_);
            y = value_type((y_ > a_) ? a_ : y_);
            k = value_type((k_ > a_) ? a_ : k_);
            a = value_type(a_);

            return *this;
        }

        //--------------------------------------------------------------------
        AGG_INLINE const self_type& demultiply()
        {
            if(a == base_mask) return *this;
            if(a == 0)
            {
                c = m = y = k = 0;
                return *this;
            }
            calc_type c_ = (calc_type(c) * base_mask) / a;
            calc_type m_ = (calc_type(m) * base_mask) / a;
            calc_type y_ = (calc_type(y) * base_mask) / a;
            calc_type k_ = (calc_type(k) * base_mask) / a;
            c = value_type((c_ > calc_type(base_mask)) ? calc_type(base_mask) : c_);
            m = value_type((m_ > calc_type(base_mask)) ? calc_type(base_mask) : m_);
            y = value_type((y_ > calc_type(base_mask)) ? calc_type(base_mask) : y_);
            k = value_type((k_ > calc_type(base_mask)) ? calc_type(base_mask) : k_);
            return *this;
        }

        //--------------------------------------------------------------------
        AGG_INLINE self_type gradient(const self_type& cl, double l) const
        {
            self_type ret;
            calc_type ik = uround(l * base_scale);
            ret.c = value_type(calc_type(c) + (((calc_type(cl.c) - c) * ik) >> base_shift));
            ret.m = value_type(calc_type(m) + (((calc_type(cl.m) - m) * ik) >> base_shift));
            ret.y = value_type(calc_type(y) + (((calc_type(cl.y) - y) * ik) >> base_shift));
            ret.k = value_type(calc_type(k) + (((calc_type(cl.k) - k) * ik) >> base_shift));
            ret.a = value_type(calc_type(a) + (((calc_type(cl.a) - a) * ik) >> base_shift));
            return ret;
        }

        //--------------------------------------------------------------------
        AGG_INLINE void add(const self_type& cl, unsigned cover)
        {
            calc_type cc, cm, cy, ck, ca;
            if(cover == cover_mask)
            {
                if(cl.a == base_mask) 
                {
                    *this = cl;
                }
                else
                {
                    cc = c + cl.c; c = (cc > calc_type(base_mask)) ? calc_type(base_mask) : cc;
                    cm = m + cl.m; m = (cm > calc_type(base_mask)) ? calc_type(base_mask) : cm;
                    cy = y + cl.y; y = (cy > calc_type(base_mask)) ? calc_type(base_mask) : cy;
                    ck = k + cl.k; k = (ck > calc_type(base_mask)) ? calc_type(base_mask) : ck;
                    ca = a + cl.a; a = (ca > calc_type(base_mask)) ? calc_type(base_mask) : ca;
                }
            }
            else
            {
                cc = c + ((cl.c * cover + cover_mask / 2) >> cover_shift);
                cm = m + ((cl.m * cover + cover_mask / 2) >> cover_shift);
                cy = y + ((cl.y * cover + cover_mask / 2) >> cover_shift);
                ck = k + ((cl.k * cover + cover_mask / 2) >> cover_shift);
                ca = a + ((cl.a * cover + cover_mask / 2) >> cover_shift);

                c = (cc > calc_type(base_mask)) ? calc_type(base_mask) : cc;
                m = (cm > calc_type(base_mask)) ? calc_type(base_mask) : cm;
                y = (cy > calc_type(base_mask)) ? calc_type(base_mask) : cy;
                k = (ck > calc_type(base_mask)) ? calc_type(base_mask) : ck;
                a = (ca > calc_type(base_mask)) ? calc_type(base_mask) : ca;
            }
        }

        //--------------------------------------------------------------------
        template<class GammaLUT>
        AGG_INLINE void apply_gamma_dir(const GammaLUT& gamma)
        {
            c = gamma.dir(c);
            m = gamma.dir(m);
            y = gamma.dir(y);
            k = gamma.dir(k);
        }

        //--------------------------------------------------------------------
        template<class GammaLUT>
        AGG_INLINE void apply_gamma_inv(const GammaLUT& gamma)
        {
            c = gamma.inv(c);
            m = gamma.inv(m);
            y = gamma.inv(y);
            k = gamma.inv(k);
        }

        //--------------------------------------------------------------------
        static self_type no_color() { return self_type(0,0,0,0); }

        static const self_type black;
        static const self_type white;
        static const self_type cyan;
        static const self_type magenta;
        static const self_type yellow;
    };
    
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
}
#endif
