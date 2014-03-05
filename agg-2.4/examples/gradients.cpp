#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "agg_rendering_buffer.h"
#include "agg_rasterizer_scanline_aa.h"
#include "agg_scanline_u.h"
#include "agg_scanline_p.h"
#include "agg_conv_transform.h"
#include "agg_color_rgba.h"
#include "agg_color_gray.h"
#include "agg_span_allocator.h"
#include "agg_span_gradient.h"
#include "agg_span_interpolator_linear.h"
#include "agg_renderer_scanline.h"
// #include "agg_path.h"
#include "ctrl/agg_rbox_ctrl.h"
#include "ctrl/agg_spline_ctrl.h"
#include "ctrl/agg_gamma_ctrl.h"
#include "platform/agg_platform_support.h"


//#define AGG_GRAY16
#define AGG_BGR24
//#define AGG_RGB24
//#define AGG_RGB_AAA
//#define AGG_RGBA32 
//#define AGG_ARGB32 
//#define AGG_ABGR32
//#define AGG_RGB565
//#define AGG_RGB555
#include "pixel_formats.h"

enum flip_y_e { flip_y = true };

const double center_x = 350;
const double center_y = 280;


class gradient_polymorphic_wrapper_base
{
public:
    virtual int calculate(int x, int y, int) const = 0;
};

namespace agg {

    /** Plugin for gradient_polymorphic_wrapper<...> which generates a
     * radial grandient as in PDF or Postscript documents.
     * The gradient is defined between two circles which do *not*
     * have to be concenctrical.
     *
     * The series of circles may be defined an a horizontal line (x0..y1)
     * with radiues (r0...r1). For a given x we nay calculate y = f(x,t) as.
     *
     * y = sqrt( ( r0 + ( r1 - r0 ) * t ) ^ 2 - ( x -  ( x0 + ( x1 - x0 ) * t ) ) ^ 2 )
     *
     * The gradient has to solve the problem t = f(x,y) for which there are two solutions.
     */

    class gradient_biradial
    {
    public:
        gradient_biradial()
            :m_r0(0),
             m_r1(0),
             m_x0(0),
             m_x1(0)
        {
        }

        ~gradient_biradial()
        {
        }


        virtual int calculate(int x,int y,int d) const
        {
            int r0 = m_r0;
            int r1 = m_r1;
            int x0 = m_x0;
            int x1 = m_x1;

            double t = f1(x,y,x0,x1,r0,r1);

            if(t!=t || t< 0 || t>1)
                t = f0(x,y,x0,x1,r0,r1);

            if(t==t && t>1)
                return d;
            if(t==t && t<0)
                return 0;
            return (t!=t || t<0 /*|| t>1*/) ? 0 : t*d;
        }

        void set_radius(int r0,int r1)
        {
            m_r0 = r0 << gradient_subpixel_shift;
            m_r1 = r1 << gradient_subpixel_shift;
        }
        void set_center(int x0,int x1)
        {
            m_x0 = x0 << gradient_subpixel_shift;
            m_x1 = x1 << gradient_subpixel_shift;
        }

    private:
        double f0(int x,int y,int x0,int x1,int r0,int r1) const
        {
            return -(::sqrt((-::pow(x1, 2) + 2 * x0 * x1 - ::pow(x0, 2) + ::pow(r1, 2) - 2 * r0 * r1 +::pow(r0, 2)) * ::pow(y, 2) + ::pow(r0, 2) * ::pow(x1, 2)
                         + ((2 * r0 * r1 - 2 * ::pow(r0, 2)) * x - 2 * r0 * r1 * x0) * x1 + ::pow(r1, 2) * ::pow(x0, 2) + (2 * r0 * r1 - 2 * ::pow(r1, 2))
                         * x * x0 + (::pow(r1, 2) - 2 * r0 * r1 + ::pow(r0, 2)) *
                           ::pow(x, 2)) + (x0 - x) * x1 - ::pow(x0, 2) + x * x0 - r0 * r1 + ::pow(r0, 2))/(::pow(x1, 2) - 2 * x0 * x1 + ::pow(x0, 2) - ::pow(r1, 2) + 2 * r0 * r1 - ::pow(r0, 2));
        }

        double f1(int x,int y,int x0,int x1,int r0,int r1) const
        {
            return (sqrt(( -::pow(x1, 2) + 2 * x0 * x1 - ::pow(x0, 2) + ::pow(r1, 2) - 2*r0*r1+::pow(r0, 2)) * ::pow(y, 2) + ::pow(r0, 2) * ::pow(x1, 2)
                             + ((2*r0*r1-2*::pow(r0, 2) ) * x - 2 * r0 * r1 * x0) * x1 + ::pow(r1, 2)
                             * ::pow(x0, 2) + ( 2 * r0 * r1 - 2 * ::pow(r1, 2) ) * x * x0
                             + ( ::pow(r1, 2) - 2 * r0 * r1 + ::pow(r0, 2)) * ::pow(x, 2)) + (x-x0) * x1 + ::pow(x0, 2) - x * x0 + r0 * r1 - ::pow(r0, 2) )
                / ( ::pow(x1, 2) - 2 * x0 * x1 + ::pow(x0, 2) - ::pow(r1, 2) + 2 * r0 * r1 - ::pow(r0, 2));
        }
        int m_r0;
        int m_r1;
        int m_x0;
        int m_x1;
    };
};

template<class GradientF> 
class gradient_polymorphic_wrapper : public gradient_polymorphic_wrapper_base
{
public:
    gradient_polymorphic_wrapper() : m_adaptor(m_gradient) {}

    virtual int calculate(int x, int y, int d) const
    {
        return m_adaptor.calculate(x, y, d);
    }
    GradientF m_gradient;
    agg::gradient_reflect_adaptor<GradientF> m_adaptor;
};



struct color_function_profile
{
    color_function_profile() {}
    color_function_profile(const color_type* colors, const agg::int8u* profile) :
        m_colors(colors), m_profile(profile) {}

    static unsigned size() { return 256; }
    const color_type& operator [] (unsigned v) const
    { 
        return m_colors[m_profile[v]]; 
    }

    const color_type* m_colors;
    const agg::int8u* m_profile;
};


class the_application : public agg::platform_support
{
    agg::gamma_ctrl<agg::rgba8>  m_profile;
    agg::spline_ctrl<agg::rgba8> m_spline_r;
    agg::spline_ctrl<agg::rgba8> m_spline_g;
    agg::spline_ctrl<agg::rgba8> m_spline_b;
    agg::spline_ctrl<agg::rgba8> m_spline_a;
    agg::rbox_ctrl<agg::rgba8>   rbox;

    double m_pdx;
    double m_pdy;
    double m_center_x;
    double m_center_y;
    double m_scale;
    double m_prev_scale;
    double m_angle;
    double m_prev_angle;
    double m_scale_x;
    double m_prev_scale_x;
    double m_scale_y;
    double m_prev_scale_y;
    bool m_mouse_move;

public:
    virtual ~the_application()
    {
        FILE* fd = fopen(full_file_name("settings.dat"), "w");
        fprintf(fd, "%f\n", m_center_x);
        fprintf(fd, "%f\n", m_center_y);
        fprintf(fd, "%f\n", m_scale);
        fprintf(fd, "%f\n", m_angle);
        fprintf(fd, "%f\n", m_spline_r.x(0));
        fprintf(fd, "%f\n", m_spline_r.y(0));
        fprintf(fd, "%f\n", m_spline_r.x(1));
        fprintf(fd, "%f\n", m_spline_r.y(1));
        fprintf(fd, "%f\n", m_spline_r.x(2));
        fprintf(fd, "%f\n", m_spline_r.y(2));
        fprintf(fd, "%f\n", m_spline_r.x(3));
        fprintf(fd, "%f\n", m_spline_r.y(3));
        fprintf(fd, "%f\n", m_spline_r.x(4));
        fprintf(fd, "%f\n", m_spline_r.y(4));
        fprintf(fd, "%f\n", m_spline_r.x(5));
        fprintf(fd, "%f\n", m_spline_r.y(5));
        fprintf(fd, "%f\n", m_spline_g.x(0));
        fprintf(fd, "%f\n", m_spline_g.y(0));
        fprintf(fd, "%f\n", m_spline_g.x(1));
        fprintf(fd, "%f\n", m_spline_g.y(1));
        fprintf(fd, "%f\n", m_spline_g.x(2));
        fprintf(fd, "%f\n", m_spline_g.y(2));
        fprintf(fd, "%f\n", m_spline_g.x(3));
        fprintf(fd, "%f\n", m_spline_g.y(3));
        fprintf(fd, "%f\n", m_spline_g.x(4));
        fprintf(fd, "%f\n", m_spline_g.y(4));
        fprintf(fd, "%f\n", m_spline_g.x(5));
        fprintf(fd, "%f\n", m_spline_g.y(5));
        fprintf(fd, "%f\n", m_spline_b.x(0));
        fprintf(fd, "%f\n", m_spline_b.y(0));
        fprintf(fd, "%f\n", m_spline_b.x(1));
        fprintf(fd, "%f\n", m_spline_b.y(1));
        fprintf(fd, "%f\n", m_spline_b.x(2));
        fprintf(fd, "%f\n", m_spline_b.y(2));
        fprintf(fd, "%f\n", m_spline_b.x(3));
        fprintf(fd, "%f\n", m_spline_b.y(3));
        fprintf(fd, "%f\n", m_spline_b.x(4));
        fprintf(fd, "%f\n", m_spline_b.y(4));
        fprintf(fd, "%f\n", m_spline_b.x(5));
        fprintf(fd, "%f\n", m_spline_b.y(5));
        fprintf(fd, "%f\n", m_spline_a.x(0));
        fprintf(fd, "%f\n", m_spline_a.y(0));
        fprintf(fd, "%f\n", m_spline_a.x(1));
        fprintf(fd, "%f\n", m_spline_a.y(1));
        fprintf(fd, "%f\n", m_spline_a.x(2));
        fprintf(fd, "%f\n", m_spline_a.y(2));
        fprintf(fd, "%f\n", m_spline_a.x(3));
        fprintf(fd, "%f\n", m_spline_a.y(3));
        fprintf(fd, "%f\n", m_spline_a.x(4));
        fprintf(fd, "%f\n", m_spline_a.y(4));
        fprintf(fd, "%f\n", m_spline_a.x(5));
        fprintf(fd, "%f\n", m_spline_a.y(5));
        double x1,y1,x2,y2;
        m_profile.values(&x1, &y1, &x2, &y2);
        fprintf(fd, "%f\n", x1);
        fprintf(fd, "%f\n", y1);
        fprintf(fd, "%f\n", x2);
        fprintf(fd, "%f\n", y2);
        fclose(fd);
    }

    the_application(agg::pix_format_e format, bool flip_y) :
        agg::platform_support(format, flip_y),
        m_profile(10.0, 10.0, 200.0, 170.0-5.0,    !flip_y),
        m_spline_r(210, 10,     210+250, 5+40,  6, !flip_y),
        m_spline_g(210, 10+40,  210+250, 5+80,  6, !flip_y),
        m_spline_b(210, 10+80,  210+250, 5+120, 6, !flip_y),
        m_spline_a(210, 10+120, 210+250, 5+160, 6, !flip_y),
        rbox(10.0, 180.0, 200.0, 320.0, !flip_y),

        m_pdx(0.0),
        m_pdy(0.0),
        m_center_x(center_x),
        m_center_y(center_y),
        m_scale(1.0),
        m_prev_scale(1.0),
        m_angle(0.0),
        m_prev_angle(0.0),
        m_scale_x(1.0),
        m_prev_scale_x(1.0),
        m_scale_y(1.0),
        m_prev_scale_y(1.0),
        m_mouse_move(false)
    {
        add_ctrl(m_profile);
        add_ctrl(m_spline_r);
        add_ctrl(m_spline_g);
        add_ctrl(m_spline_b);
        add_ctrl(m_spline_a);
        add_ctrl(rbox);

        m_profile.border_width(2.0, 2.0);

        m_spline_r.background_color(agg::rgba(1.0, 0.8, 0.8));
        m_spline_g.background_color(agg::rgba(0.8, 1.0, 0.8));
        m_spline_b.background_color(agg::rgba(0.8, 0.8, 1.0));
        m_spline_a.background_color(agg::rgba(1.0, 1.0, 1.0));

        m_spline_r.border_width(1.0, 2.0);
        m_spline_g.border_width(1.0, 2.0);
        m_spline_b.border_width(1.0, 2.0);
        m_spline_a.border_width(1.0, 2.0);
        rbox.border_width(2.0, 2.0);

        m_spline_r.point(0, 0.0,     1.0);
        m_spline_r.point(1, 1.0/5.0, 1.0 - 1.0/5.0);
        m_spline_r.point(2, 2.0/5.0, 1.0 - 2.0/5.0);
        m_spline_r.point(3, 3.0/5.0, 1.0 - 3.0/5.0);
        m_spline_r.point(4, 4.0/5.0, 1.0 - 4.0/5.0);
        m_spline_r.point(5, 1.0,     0.0);
        m_spline_r.update_spline();

        m_spline_g.point(0, 0.0,     1.0);
        m_spline_g.point(1, 1.0/5.0, 1.0 - 1.0/5.0);
        m_spline_g.point(2, 2.0/5.0, 1.0 - 2.0/5.0);
        m_spline_g.point(3, 3.0/5.0, 1.0 - 3.0/5.0);
        m_spline_g.point(4, 4.0/5.0, 1.0 - 4.0/5.0);
        m_spline_g.point(5, 1.0,     0.0);
        m_spline_g.update_spline();

        m_spline_b.point(0, 0.0,     1.0);
        m_spline_b.point(1, 1.0/5.0, 1.0 - 1.0/5.0);
        m_spline_b.point(2, 2.0/5.0, 1.0 - 2.0/5.0);
        m_spline_b.point(3, 3.0/5.0, 1.0 - 3.0/5.0);
        m_spline_b.point(4, 4.0/5.0, 1.0 - 4.0/5.0);
        m_spline_b.point(5, 1.0,     0.0);
        m_spline_b.update_spline();

        m_spline_a.point(0, 0.0,     1.0);
        m_spline_a.point(1, 1.0/5.0, 1.0);
        m_spline_a.point(2, 2.0/5.0, 1.0);
        m_spline_a.point(3, 3.0/5.0, 1.0);
        m_spline_a.point(4, 4.0/5.0, 1.0);
        m_spline_a.point(5, 1.0,     1.0);
        m_spline_a.update_spline();

        rbox.add_item("Bicircular");
        rbox.add_item("Circular");
        rbox.add_item("Diamond");
        rbox.add_item("Linear");
        rbox.add_item("XY");
        rbox.add_item("sqrt(XY)");
        rbox.add_item("Conic");
        rbox.cur_item(0);

        FILE* fd = fopen(full_file_name("settings.dat"), "r");

        if(fd)
        {
            float x;
            float y;
            float x2;
            float y2;
            float t;

            fscanf(fd, "%f\n", &t); m_center_x = t;
            fscanf(fd, "%f\n", &t); m_center_y = t;
            fscanf(fd, "%f\n", &t); m_scale = t;
            fscanf(fd, "%f\n", &t); m_angle = t;
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_r.point(0, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_r.point(1, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_r.point(2, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_r.point(3, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_r.point(4, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_r.point(5, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_g.point(0, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_g.point(1, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_g.point(2, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_g.point(3, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_g.point(4, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_g.point(5, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_b.point(0, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_b.point(1, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_b.point(2, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_b.point(3, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_b.point(4, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_b.point(5, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_a.point(0, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_a.point(1, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_a.point(2, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_a.point(3, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_a.point(4, x, y);
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y); m_spline_a.point(5, x, y);
            m_spline_r.update_spline();
            m_spline_g.update_spline();
            m_spline_b.update_spline();
            m_spline_a.update_spline();
            fscanf(fd, "%f\n", &x);
            fscanf(fd, "%f\n", &y);
            fscanf(fd, "%f\n", &x2);
            fscanf(fd, "%f\n", &y2);
            m_profile.values(x, y, x2, y2);
            fclose(fd);
        }

    }


    virtual void on_draw()
    {
        agg::rasterizer_scanline_aa<> ras;

        typedef agg::renderer_base<pixfmt> renderer_base;
        agg::scanline_u8 sl;


        pixfmt pixf(rbuf_window());
        renderer_base rb(pixf);
        rb.clear(agg::rgba(0, 0, 0));
      
        m_profile.text_size(8.0);

        agg::render_ctrl(ras, sl, rb, m_profile);
        agg::render_ctrl(ras, sl, rb, m_spline_r);
        agg::render_ctrl(ras, sl, rb, m_spline_g);
        agg::render_ctrl(ras, sl, rb, m_spline_b);
        agg::render_ctrl(ras, sl, rb, m_spline_a);
        agg::render_ctrl(ras, sl, rb, rbox);

        double ini_scale = 1.0;

        agg::trans_affine mtx1;

        mtx1 *= agg::trans_affine_scaling(ini_scale, ini_scale);
        mtx1 *= agg::trans_affine_rotation(agg::deg2rad(0.0));
        mtx1 *= agg::trans_affine_translation(center_x, center_y);
        mtx1 *= trans_affine_resizing();
#if 1
        agg::ellipse e1;
        e1.init(0.0, 0.0, 110.0, 110.0, 64);

#else
        agg::path_storage ps;

        ps.move_to(-110,-110);
        ps.line_to(-110,110);
        ps.line_to(110,110);
        ps.line_to(110,-110);
        ps.line_to(-110,-110);
        ps.close_polygon();
#endif

        agg::trans_affine mtx_g1;

        mtx_g1 *= agg::trans_affine_scaling(ini_scale, ini_scale);
        mtx_g1 *= agg::trans_affine_scaling(m_scale, m_scale);
        mtx_g1 *= agg::trans_affine_scaling(m_scale_x, m_scale_y);
        mtx_g1 *= agg::trans_affine_rotation( m_angle );
        mtx_g1 *= agg::trans_affine_translation(m_center_x, m_center_y);
        // mtx_g1 *= agg::trans_affine_translation(center_x-110, center_y-110);
        mtx_g1 *= trans_affine_resizing();


        mtx_g1.invert();


        color_type color_profile[256]; // color_type is defined in pixel_formats.h
        int i;
        for(i = 0; i < 256; i++)
        {
            color_profile[i] = color_type(agg::rgba(m_spline_r.spline()[i], 
                                                    m_spline_g.spline()[i],
                                                    m_spline_b.spline()[i],
                                                    m_spline_a.spline()[i]));
        }
#if 1
        agg::conv_transform<agg::ellipse, agg::trans_affine> t1(e1, mtx1);
#else
        agg::conv_transform<agg::path_storage, agg::trans_affine> t1(ps, mtx1);
#endif
        gradient_polymorphic_wrapper<agg::gradient_biradial>     gr_bicircle;
        gradient_polymorphic_wrapper<agg::gradient_radial>       gr_circle;
        gradient_polymorphic_wrapper<agg::gradient_diamond>      gr_diamond;
        gradient_polymorphic_wrapper<agg::gradient_x>            gr_x;
        gradient_polymorphic_wrapper<agg::gradient_xy>           gr_xy;
        gradient_polymorphic_wrapper<agg::gradient_sqrt_xy>      gr_sqrt_xy;
        gradient_polymorphic_wrapper<agg::gradient_conic>        gr_conic;

        gr_bicircle.m_gradient.set_radius(15,75);
        gr_bicircle.m_gradient.set_center(0,120);

        gradient_polymorphic_wrapper_base* gr_ptr = &gr_bicircle;

        // gr_circle.m_gradient.init(150, 80, 80);

        switch(rbox.cur_item())
        {
            case 1: gr_ptr = &gr_circle;   break;
            case 2: gr_ptr = &gr_diamond; break;
            case 3: gr_ptr = &gr_x;       break;
            case 4: gr_ptr = &gr_xy;      break;
            case 5: gr_ptr = &gr_sqrt_xy; break;
            case 6: gr_ptr = &gr_conic;   break;
        }

        typedef agg::span_interpolator_linear<> interpolator_type;
        typedef agg::span_gradient<color_type, 
                                   interpolator_type,
                                   gradient_polymorphic_wrapper_base,
                                   color_function_profile> gradient_span_gen;
        typedef agg::span_allocator<gradient_span_gen::color_type> gradient_span_alloc;

        gradient_span_alloc    span_alloc;
        color_function_profile colors(color_profile, m_profile.gamma());
        interpolator_type      inter(mtx_g1);

        gradient_span_gen      span_gen(inter, *gr_ptr, colors, 0, 150);

        ras.add_path(t1);

        agg::render_scanlines_aa(ras, sl, rb, span_alloc, span_gen);
    }


    virtual void on_mouse_move(int x, int y, unsigned flags)
    {
        if(m_mouse_move)
        {
            double x2 = x;
            double y2 = y;
            trans_affine_resizing().inverse_transform(&x2, &y2);

            if(flags & agg::kbd_ctrl)
            {
                double dx = x2 - m_center_x;
                double dy = y2 - m_center_y;
                m_scale_x = m_prev_scale_x * dx / m_pdx;
                m_scale_y = m_prev_scale_y * dy / m_pdy;
                force_redraw();
            }
            else
            {
                if(flags & agg::mouse_left)
                {
                    m_center_x = x2 + m_pdx;
                    m_center_y = y2 + m_pdy;
                    force_redraw();
                }

                if(flags & agg::mouse_right)
                {
                    double dx = x2 - m_center_x;
                    double dy = y2 - m_center_y;
                    m_scale = m_prev_scale * 
                              sqrt(dx * dx + dy * dy) / 
                              sqrt(m_pdx * m_pdx + m_pdy * m_pdy);

                    m_angle = m_prev_angle + atan2(dy, dx) - atan2(m_pdy, m_pdx);
                    force_redraw();
                }
            }
        }
    }


    virtual void on_mouse_button_down(int x, int y, unsigned flags)
    {
        m_mouse_move = true;
        double x2 = x;
        double y2 = y;
        trans_affine_resizing().inverse_transform(&x2, &y2);

        m_pdx = m_center_x - x2;
        m_pdy = m_center_y - y2;
        m_prev_scale = m_scale;
        m_prev_angle = m_angle + agg::pi;
        m_prev_scale_x = m_scale_x;
        m_prev_scale_y = m_scale_y;
        force_redraw();
    }


    virtual void on_mouse_button_up(int x, int y, unsigned flags)
    {
        m_mouse_move = false;
    }


    virtual void on_key(int x, int y, unsigned key, unsigned flags) 
    {
        if(key == agg::key_f1)
        {
            FILE* fd = fopen(full_file_name("colors.dat"), "w");
            int i;
            for(i = 0; i < 256; i++)
            {
                color_type c = agg::rgba(m_spline_r.spline()[i], 
                                         m_spline_g.spline()[i],
                                         m_spline_b.spline()[i],
                                         m_spline_a.spline()[i]);
                fprintf(fd, "    %3d, %3d, %3d, %3d,\n", c.r, c.g, c.b, c.a);
            }
            fclose(fd);

            fd = fopen(full_file_name("profile.dat"), "w");
            for(i = 0; i < 256; i++)
            {
                fprintf(fd, "%3d, ", unsigned(m_profile.gamma()[i]));
                if((i & 0xF) == 0xF) fprintf(fd, "\n");
            }
            fclose(fd);
        }


    }

};



int agg_main(int argc, char* argv[])
{
    //#ifdef _WIN32
    //    FILE* fd = fopen("stdout.txt", "w"); fclose(fd);
    //#endif
    //AGG_WATCHDOGGY(wd1, false);

    the_application app(pix_format, flip_y);
    app.caption("AGG gradients with Mach bands compensation");

    if(app.init(512, 400, agg::window_resize | agg::window_hw_buffer))
    {
        return app.run();
    }
    return 1;
}


