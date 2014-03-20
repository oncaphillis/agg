//----------------------------------------------------------------------------
// Anti-Grain Geometry - Version 2.4
// Copyright (C) 2002-2005 Maxim Shemanarev (http://www.antigrain.com)
// Copyright (C) 2014 Sebastian Kloska (oncaphillis@snafu.de)
//
// Permission to copy, use, modify, sell and distribute this software 
// is granted provided this copyright notice appears in all copies. 
// This software is provided "as is" without express or implied
// warranty, and with no claim as to its suitability for any purpose.
//
//----------------------------------------------------------------------------
// Contact: mcseem@antigrain.com
//          mcseemagg@yahoo.com
//          http://www.antigrain.com
//----------------------------------------------------------------------------

#ifndef AGG_SPAN_GRADIENT_INCLUDED
#define AGG_SPAN_GRADIENT_INCLUDED

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "agg_basics.h"
#include "agg_math.h"
#include "agg_array.h"


namespace agg
{
    // A bit if SFINAE magic helps us to make the subpixel shift value
    // variable. This value was previously defined as the enum
    // gradient_subpixel_shift and determines the size of the precision
    // for AGGs fixed-point arithmetic. The previous value was 4 which
    // allows for 16 different values after the dot. In some cases this might
    // result in jaggy gradients.
    // The SFINAE struct allows the value of subpixel_shift to be determined
    // via the concrete GradientF template argument of span_gradient<...>
    // If the class does not provide its own constant the old constant value
    // will still be used. Otherwise we "inherit" the value from the gradient function.

    template<class GradientF>
    class gradient_props
    {
    private:
        typedef char Yes;
        typedef char No[2];

        template<typename X,const X *> struct has;
        template <typename Y> static Yes & Test(has<int,&Y::subpixel_shift>*);
        template <typename Y> static  No & Test(...);
        template<class X,bool B>
        struct V
        {
            static const int value = 4;
        };

        template<class X>
        struct V<X,true>
        {
            static const int value=X::subpixel_shift;
        };
    public:
        static const int subpixel_shift = V<GradientF,sizeof(Test<GradientF>(0))==sizeof(Yes)>::value;
        static const int subpixel_scale = 1 << subpixel_shift;
        static const int subpixel_mask  = subpixel_scale - 1;
    };

    //==========================================================span_gradient
    template<class ColorT,
             class Interpolator,
             class GradientF, 
             class ColorF>
    class span_gradient
    {
    public:
        typedef Interpolator interpolator_type;
        typedef ColorT color_type;

        enum downscale_shift_e
        {
            downscale_shift = interpolator_type::subpixel_shift - gradient_props<GradientF>::subpixel_shift
        };

        //--------------------------------------------------------------------
        span_gradient() {}

        //--------------------------------------------------------------------
        span_gradient(interpolator_type& inter,
                      const GradientF& gradient_function,
                      const ColorF& color_function,
                      double d1, double d2) : 
            m_interpolator(&inter),
            m_gradient_function(&gradient_function),
            m_color_function(&color_function),
            m_d1(iround(d1 * gradient_props<GradientF>::subpixel_scale)),
            m_d2(iround(d2 * gradient_props<GradientF>::subpixel_scale))
        {
        }

        //--------------------------------------------------------------------
        interpolator_type& interpolator() { return *m_interpolator; }
        const GradientF& gradient_function() const { return *m_gradient_function; }
        const ColorF& color_function() const { return *m_color_function; }
        double d1() const { return double(m_d1) / gradient_props<GradientF>::subpixel_scale; }
        double d2() const { return double(m_d2) / gradient_props<GradientF>::subpixel_scale; }

        //--------------------------------------------------------------------
        void interpolator(interpolator_type& i) { m_interpolator = &i; }
        void gradient_function(const GradientF& gf) { m_gradient_function = &gf; }
        void color_function(const ColorF& cf) { m_color_function = &cf; }
        void d1(double v) { m_d1 = iround(v * gradient_props<GradientF>::subpixel_scale); }
        void d2(double v) { m_d2 = iround(v * gradient_props<GradientF>::subpixel_scale); }

        //--------------------------------------------------------------------
        void prepare() {}

        //--------------------------------------------------------------------
        void generate(color_type* span, int x, int y, unsigned len)
        {   
            int dd = m_d2 - m_d1;
            if(dd < 1) dd = 1;
            m_interpolator->begin(x+0.5, y+0.5, len);
            do
            {
                m_interpolator->coordinates(&x, &y);
                int d = m_gradient_function->calculate(x >> downscale_shift, 
                                                       y >> downscale_shift, m_d2);
                d = ((d - m_d1) * (int)m_color_function->size()) / dd;
                if(d < 0) d = 0;
                if(d >= (int)m_color_function->size()) d = m_color_function->size() - 1;
                *span++ = (*m_color_function)[d];
                ++(*m_interpolator);
            }
            while(--len);
        }

    private:
        interpolator_type* m_interpolator;
        const GradientF*   m_gradient_function;
        const ColorF*      m_color_function;
        int                m_d1;
        int                m_d2;
    };




    //=====================================================gradient_linear_color
    template<class ColorT> 
    struct gradient_linear_color
    {
        typedef ColorT color_type;

        gradient_linear_color() {}
        gradient_linear_color(const color_type& c1, const color_type& c2, 
                              unsigned size = 256) :
            m_c1(c1), m_c2(c2), m_size(size)
				// VFALCO 4/28/09
				,m_mult(1/(double(size)-1))
				// VFALCO
			{}

        unsigned size() const { return m_size; }
        color_type operator [] (unsigned v) const 
        {
			// VFALCO 4/28/09 
            //return m_c1.gradient(m_c2, double(v) / double(m_size - 1));
            return m_c1.gradient(m_c2, double(v) * m_mult );
			// VFALCO
        }

        void colors(const color_type& c1, const color_type& c2, unsigned size = 256)
        {
            m_c1 = c1;
            m_c2 = c2;
            m_size = size;
			// VFALCO 4/28/09
			m_mult=1/(double(size)-1);
			// VFALCO
        }

        color_type m_c1;
        color_type m_c2;
        unsigned m_size;
		// VFALCO 4/28/09
		double m_mult;
		// VFALCO
    };

    //==========================================================gradient_circle
    class gradient_circle
    {
        // Actually the same as radial. Just for compatibility
    public:
        static AGG_INLINE int calculate(int x, int y, int)
        {
            return int(fast_sqrt(x*x + y*y));
        }
    };


    //==========================================================gradient_radial
    class gradient_radial
    {
    public:
        static AGG_INLINE int calculate(int x, int y, int)
        {
            return int(fast_sqrt(x*x + y*y));
        }
    };

    //========================================================gradient_radial_d
    class gradient_radial_d
    {
    public:
        static AGG_INLINE int calculate(int x, int y, int)
        {
            return uround(sqrt(double(x)*double(x) + double(y)*double(y)));
        }
    };

    //====================================================gradient_radial_focus
    template<class T=void>
    class basic_gradient_radial_focus
    {
    public:
        
        //---------------------------------------------------------------------
        basic_gradient_radial_focus() :
            m_r(100 * gradient_props<T>::subpixel_scale),
            m_fx(0), 
            m_fy(0)
        {
            update_values();
        }

        //---------------------------------------------------------------------
        basic_gradient_radial_focus(double r, double fx, double fy) :
            m_r (iround(r  * gradient_props<T>::subpixel_scale)),
            m_fx(iround(fx * gradient_props<T>::subpixel_scale)),
            m_fy(iround(fy * gradient_props<T>::subpixel_scale))
        {
            update_values();
        }

        //---------------------------------------------------------------------
        void init(double r, double fx, double fy)
        {
            m_r  = iround(r  * gradient_props<T>::subpixel_scale);
            m_fx = iround(fx * gradient_props<T>::subpixel_scale);
            m_fy = iround(fy * gradient_props<T>::subpixel_scale);
            update_values();
        }

        //---------------------------------------------------------------------
        double radius()  const { return double(m_r)  / gradient_props<T>::subpixel_scale; }
        double focus_x() const { return double(m_fx) / gradient_props<T>::subpixel_scale; }
        double focus_y() const { return double(m_fy) / gradient_props<T>::subpixel_scale; }

        //---------------------------------------------------------------------
        int calculate(int x, int y, int) const
        {
            double dx = x - m_fx;
            double dy = y - m_fy;
            double d2 = dx * m_fy - dy * m_fx;
            double d3 = m_r2 * (dx * dx + dy * dy) - d2 * d2;
            return iround((dx * m_fx + dy * m_fy + sqrt(fabs(d3))) * m_mul);
        }

    private:
        //---------------------------------------------------------------------
        void update_values()
        {
            // Calculate the invariant values. In case the focal center
            // lies exactly on the gradient circle the divisor degenerates
            // into zero. In this case we just move the focal center by
            // one subpixel unit possibly in the direction to the origin (0,0)
            // and calculate the values again.
            //-------------------------
            m_r2  = double(m_r)  * double(m_r);
            m_fx2 = double(m_fx) * double(m_fx);
            m_fy2 = double(m_fy) * double(m_fy);
            double d = (m_r2 - (m_fx2 + m_fy2));
            if(d == 0)
            {
                if(m_fx) { if(m_fx < 0) ++m_fx; else --m_fx; }
                if(m_fy) { if(m_fy < 0) ++m_fy; else --m_fy; }
                m_fx2 = double(m_fx) * double(m_fx);
                m_fy2 = double(m_fy) * double(m_fy);
                d = (m_r2 - (m_fx2 + m_fy2));
            }
            m_mul = m_r / d;
        }

        int    m_r;
        int    m_fx;
        int    m_fy;
        double m_r2;
        double m_fx2;
        double m_fy2;
        double m_mul;
    };
    typedef basic_gradient_radial_focus<> gradient_radial_focus;

    //==============================================================gradient_biradial

    /** A gradient which generates a
     * radial grandient as in PDF or Postscript documents.
     * The gradient is defined between two circles which do *not*
     * have to be concenctrical.
     *
     * The series of circles may be defined an a horizontal line (x0..y1)
     * with radiues (r0...r1). For a given x we nay calculate y = f(x,t) as.
     *
     * y = sqrt( ( r0 + ( r1 - r0 ) * t ) ^ 2 - ( x -  ( x0 + ( x1 - x0 ) * t ) ) ^ 2 )
     *
     * The gradient has to solve the problem t = f(x,y) for which there are two solutions
     * implemented in the method fl(...) and fr(...)
     *
     * In contrast to theother gradient_* classes the calculate method of this class
     * is *not* static since we need the values r0 and r1 attached to the gradient.
     */

    template<typename T=void>
    class basic_gradient_biradial
    {
    public:
        basic_gradient_biradial()
            :m_r0(0), m_r1(0), m_x0(0), m_x1(0), m_xi(0),
             m_c0(0), m_c1(0), m_c2(0), m_c3(0),
             m_c4(0), m_c5(0), m_c6(0)
        {
        }

        ~basic_gradient_biradial()
        {
        }

        int calculate(int x,int y,int d) const
        {
            double t = fl(x,y);

            if( t!=t || t> 1 )
                t = fr(x,y);
            
            // Here we have to fine tune whenever we implment a proper [extend,extend]
            // like in PDF or Postscript
            return t!=t ? 0 : t * d < m_xi ? 0 : t < 0 ? 0 : t > 1 ? d : t * d;
        }

        void set_radius(int r0,int r1)
        {
            m_r0 = r0 << gradient_props<T>::subpixel_shift;
            m_r1 = r1 << gradient_props<T>::subpixel_shift;
            calc_const();
        }

        void set_center(int x0,int x1)
        {
            m_x0 = x0 << gradient_props<T>::subpixel_shift;
            m_x1 = x1 << gradient_props<T>::subpixel_shift;
            calc_const();
        }

    private:

        /** Recalculate the constants for our t=f(x,y). There are a couple
         * of values only depending on r0..r1 and x0..x1 which we only have
         * to calculate once whenever these values change.
         */

        void calc_const()
        {
            m_c0 =  ( m_x1 * m_x1) - 2 * m_x0 * m_x1 + (m_x0 * m_x0) - (m_r1 * m_r1) + 2 * m_r0 * m_r1 - (m_r0 * m_r0);
            m_c1 = -( m_x1 * m_x1) + 2 * m_x0 * m_x1 - (m_x0 * m_x0) + (m_r1 * m_r1) - 2 * m_r0 * m_r1 + (m_r0 * m_r0);
            m_c2 = m_r0 * m_r0 * m_x1 * m_x1;
            m_c3 = 2 * m_r0 * m_r1 - 2 * m_r0 * m_r0;
            m_c4 = 2 * m_r0 * m_r1 * m_x0;
            m_c5 = (2 * m_r0 * m_r1 - 2 * (m_r1 * m_r1));
            m_c6 = (m_r1 * m_r1) - 2 * m_r0 * m_r1 + (m_r0* m_r0);

            // This is the absolut minimum for which we calculate values
            m_xi = -m_r0 / ((m_r1 - m_r0) / (m_x1 - m_x0)) + m_x0;
        }

        // The first solution for t = f(x,y)
        double fr(int x,int y) const
        {
            return -(::sqrt( double(m_c1) * double(y *y) + m_c2 + (m_c3 * x - m_c4)
                         * m_x1 + (m_r1 * m_r1) * (m_x0 * m_x0) + m_c5 * x * m_x0 + m_c6 * double(x * x))
                      + (m_x0 - x) * m_x1 - (m_x0 * m_x0) + double(x) * m_x0 - m_r0 * m_r1 + (m_r0 * m_r0))
                    / double(m_c0);
        }

        // The second solution for t = f(x,y)
        double fl(int x,int y) const
        {
            return (sqrt( double(m_c1) * double(y * y) + m_c2 + (m_c3 * x - m_c4) * m_x1
                          + (m_r1 * m_r1) * (m_x0 * m_x0) + m_c5 * x * m_x0 + m_c6 * double(x*x))
                    + double(x-m_x0) * m_x1 + (m_x0 * m_x0) - double(x) * m_x0 + m_r0 * m_r1 - (m_r0 * m_r0) )
                / double(m_c0);
        }

        // minimal and maximal x and r
        double m_r0, m_r1, m_x0, m_x1,m_xi;

        // various constants we can pre-calc for t = f(x,y)
        double m_c0, m_c1, m_c2, m_c3;
        double m_c4, m_c5, m_c6;
    };

    typedef basic_gradient_biradial<> gradient_biradial;

    //==============================================================gradient_x
    class gradient_x
    {
    public:
        static int calculate(int x, int, int) { return x; }
    };


    //==============================================================gradient_y
    class gradient_y
    {
    public:
        static int calculate(int, int y, int) { return y; }
    };

    //========================================================gradient_diamond
    class gradient_diamond
    {
    public:
        static AGG_INLINE int calculate(int x, int y, int) 
        { 
            int ax = abs(x);
            int ay = abs(y);
            return ax > ay ? ax : ay; 
        }
    };

    //=============================================================gradient_xy
    class gradient_xy
    {
    public:
        static AGG_INLINE int calculate(int x, int y, int d) 
        { 
            return abs(x) * abs(y) / d; 
        }
    };

    //========================================================gradient_sqrt_xy
    class gradient_sqrt_xy
    {
    public:
        static AGG_INLINE int calculate(int x, int y, int) 
        { 
            return fast_sqrt(abs(x) * abs(y)); 
        }
    };

    //==========================================================gradient_conic
    class gradient_conic
    {
    public:
        static AGG_INLINE int calculate(int x, int y, int d) 
        { 
            return uround(fabs(atan2(double(y), double(x))) * double(d) / pi);
        }
    };

    //=================================================gradient_repeat_adaptor
    template<class GradientF> class gradient_repeat_adaptor
    {
    public:
        gradient_repeat_adaptor(const GradientF& gradient) : 
            m_gradient(&gradient) {}

        AGG_INLINE int calculate(int x, int y, int d) const
        {
            int ret = m_gradient->calculate(x, y, d) % d;
            if(ret < 0) ret += d;
            return ret;
        }

    private:
        const GradientF* m_gradient;
    };

    //================================================gradient_reflect_adaptor
    template<class GradientF> class gradient_reflect_adaptor
    {
    public:

        static const int subpixel_shift = gradient_props<GradientF>::subpixel_shift;
        static const int subpixel_scale = gradient_props<GradientF>::subpixel_scale;

        gradient_reflect_adaptor(const GradientF& gradient) : 
            m_gradient(&gradient) {}

        AGG_INLINE int calculate(int x, int y, int d) const
        {
            int d2 = d << 1;
            int ret = m_gradient->calculate(x, y, d) % d2;
            if(ret <  0) ret += d2;
            if(ret >= d) ret  = d2 - ret;
            return ret;
        }

    private:
        const GradientF* m_gradient;
    };


}

#endif
