//----------------------------------------------------------------------------
// Anti-Grain Geometry - Version 2.4
// Copyright (C) 2002-2005 Maxim Shemanarev (http://www.antigrain.com)
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
//
// Adaptation for high precision colors has been sponsored by 
// Liberty Technology Systems, Inc., visit http://lib-sys.com
//
// Liberty Technology Systems, Inc. is the provider of
// PostScript and PDF technology for software developers.
// 
//----------------------------------------------------------------------------

#ifndef AGG_PIXFMT_RGB_INCLUDED
#define AGG_PIXFMT_RGB_INCLUDED

#include <string.h>
#include "agg_basics.h"
#include "agg_color_rgba.h"
#include "agg_rendering_buffer.h"

namespace agg
{

    //=====================================================apply_gamma_dir_rgb
    template<class ColorT, class Order, class GammaLut> class apply_gamma_dir_rgb
    {
    public:
        typedef typename ColorT::value_type value_type;

        apply_gamma_dir_rgb(const GammaLut& gamma) : m_gamma(gamma) {}

        AGG_INLINE void operator () (value_type* p)
        {
            p[Order::R] = m_gamma.dir(p[Order::R]);
            p[Order::G] = m_gamma.dir(p[Order::G]);
            p[Order::B] = m_gamma.dir(p[Order::B]);
        }

    private:
        const GammaLut& m_gamma;
    };



    //=====================================================apply_gamma_inv_rgb
    template<class ColorT, class Order, class GammaLut> class apply_gamma_inv_rgb
    {
    public:
        typedef typename ColorT::value_type value_type;

        apply_gamma_inv_rgb(const GammaLut& gamma) : m_gamma(gamma) {}

        AGG_INLINE void operator () (value_type* p)
        {
            p[Order::R] = m_gamma.inv(p[Order::R]);
            p[Order::G] = m_gamma.inv(p[Order::G]);
            p[Order::B] = m_gamma.inv(p[Order::B]);
        }

    private:
        const GammaLut& m_gamma;
    };


    //=========================================================blender_rgb
    template<class ColorT, class Order> struct blender_rgb
    {
        typedef ColorT color_type;
        typedef Order order_type;
        typedef typename color_type::value_type value_type;
        typedef typename color_type::calc_type calc_type;
        enum base_scale_e { base_shift = color_type::base_shift };

        // Blend pixels using the non-premultiplied form of Alvy-Ray Smith's
        // compositing function. Since the render buffer is opaque we skip the
        // initial premultiply and final demultiply.

        //--------------------------------------------------------------------
        static AGG_INLINE void blend_pix(value_type* p, 
                                         unsigned cr, unsigned cg, unsigned cb, 
                                         unsigned alpha, 
                                         unsigned cover)
        {
            alpha = color_type::int_mult_cover(alpha, cover);
            p[Order::R] = (value_type)(color_type::int_lerp(p[Order::R], cr, alpha));
            p[Order::G] = (value_type)(color_type::int_lerp(p[Order::G], cg, alpha));
            p[Order::B] = (value_type)(color_type::int_lerp(p[Order::B], cb, alpha));
        }
        
        //--------------------------------------------------------------------
        static AGG_INLINE void blend_pix(value_type* p, 
                                         unsigned cr, unsigned cg, unsigned cb, 
                                         unsigned alpha)
        {
            p[Order::R] = (value_type)(color_type::int_lerp(p[Order::R], cr, alpha));
            p[Order::G] = (value_type)(color_type::int_lerp(p[Order::G], cg, alpha));
            p[Order::B] = (value_type)(color_type::int_lerp(p[Order::B], cb, alpha));
        }
    };


    //======================================================blender_rgb_pre
    template<class ColorT, class Order> struct blender_rgb_pre
    {
        typedef ColorT color_type;
        typedef Order order_type;
        typedef typename color_type::value_type value_type;
        typedef typename color_type::calc_type calc_type;
        enum base_scale_e { base_shift = color_type::base_shift };

        // Blend pixels using the premultiplied form of Alvy-Ray Smith's
        // compositing function. 

        //--------------------------------------------------------------------
        static AGG_INLINE void blend_pix(value_type* p, 
                                         unsigned cr, unsigned cg, unsigned cb,
                                         unsigned alpha,
                                         unsigned cover)
        {
            cr = color_type::int_mult_cover(cr, cover);
            cg = color_type::int_mult_cover(cg, cover);
            cb = color_type::int_mult_cover(cb, cover);
            alpha = color_type::int_mult_cover(alpha, cover);
            p[Order::R] = (value_type)(color_type::int_prelerp(p[Order::R], cr, alpha));
            p[Order::G] = (value_type)(color_type::int_prelerp(p[Order::G], cg, alpha));
            p[Order::B] = (value_type)(color_type::int_prelerp(p[Order::B], cb, alpha));
        }

        //--------------------------------------------------------------------
        static AGG_INLINE void blend_pix(value_type* p, 
                                         unsigned cr, unsigned cg, unsigned cb,
                                         unsigned alpha)
        {
            p[Order::R] = (value_type)(color_type::int_prelerp(p[Order::R], cr, alpha));
            p[Order::G] = (value_type)(color_type::int_prelerp(p[Order::G], cg, alpha));
            p[Order::B] = (value_type)(color_type::int_prelerp(p[Order::B], cb, alpha));
        }

    };



    //===================================================blender_rgb_gamma
    template<class ColorT, class Order, class Gamma> class blender_rgb_gamma
    {
    public:
        typedef ColorT color_type;
        typedef Order order_type;
        typedef Gamma gamma_type;
        typedef typename color_type::value_type value_type;
        typedef typename color_type::calc_type calc_type;
        enum base_scale_e { base_shift = color_type::base_shift };

        //--------------------------------------------------------------------
        blender_rgb_gamma() : m_gamma(0) {}
        void gamma(const gamma_type& g) { m_gamma = &g; }

        //--------------------------------------------------------------------
        AGG_INLINE void blend_pix(value_type* p, 
                                  unsigned cr, unsigned cg, unsigned cb,
                                  unsigned alpha, 
                                  unsigned cover)
        {
            alpha = color_type::int_mult_cover(alpha, cover);
            calc_type r = m_gamma->dir(p[Order::R]);
            calc_type g = m_gamma->dir(p[Order::G]);
            calc_type b = m_gamma->dir(p[Order::B]);
            p[Order::R] = m_gamma->inv(color_type::int_lerp(r, m_gamma->dir(cr), alpha));
            p[Order::G] = m_gamma->inv(color_type::int_lerp(g, m_gamma->dir(cg), alpha));
            p[Order::B] = m_gamma->inv(color_type::int_lerp(b, m_gamma->dir(cb), alpha));
        }
        
        //--------------------------------------------------------------------
        AGG_INLINE void blend_pix(value_type* p, 
                                  unsigned cr, unsigned cg, unsigned cb,
                                  unsigned alpha)
        {
            calc_type r = m_gamma->dir(p[Order::R]);
            calc_type g = m_gamma->dir(p[Order::G]);
            calc_type b = m_gamma->dir(p[Order::B]);
            p[Order::R] = m_gamma->inv(color_type::int_lerp(r, m_gamma->dir(cr), alpha));
            p[Order::G] = m_gamma->inv(color_type::int_lerp(g, m_gamma->dir(cg), alpha));
            p[Order::B] = m_gamma->inv(color_type::int_lerp(b, m_gamma->dir(cb), alpha));
        }
        
    private:
        const gamma_type* m_gamma;
    };



    
    //==================================================pixfmt_alpha_blend_rgb
    template<class Blender, class RenBuf, const int Step = 3> class pixfmt_alpha_blend_rgb
    {
    public:
        typedef RenBuf   rbuf_type;
        typedef Blender  blender_type;
        typedef typename rbuf_type::row_data row_data;
        typedef typename blender_type::color_type color_type;
        typedef typename blender_type::order_type order_type;
        typedef typename color_type::value_type value_type;
        typedef typename color_type::calc_type calc_type;
        enum base_scale_e 
        {
            base_shift = color_type::base_shift,
            base_scale = color_type::base_scale,
            base_mask  = color_type::base_mask,
            pix_width  = sizeof(value_type) * Step
        };

    private:
        //--------------------------------------------------------------------
        AGG_INLINE void copy_or_blend_pix(value_type* p, 
                                          const color_type& c, 
                                          unsigned cover)
        {
            if (c.a)
            {
                if(c.a == base_mask && cover == cover_mask)
                {
                    p[order_type::R] = c.r;
                    p[order_type::G] = c.g;
                    p[order_type::B] = c.b;
                }
                else
                {
                    m_blender.blend_pix(p, c.r, c.g, c.b, c.a, cover);
                }
            }
        }

        //--------------------------------------------------------------------
        AGG_INLINE void copy_or_blend_pix(value_type* p, 
                                          const color_type& c)
        {
            if (c.a)
            {
                if(c.a == base_mask)
                {
                    p[order_type::R] = c.r;
                    p[order_type::G] = c.g;
                    p[order_type::B] = c.b;
                }
                else
                {
                    m_blender.blend_pix(p, c.r, c.g, c.b, c.a);
                }
            }
        }


    public:
        //--------------------------------------------------------------------
        explicit pixfmt_alpha_blend_rgb(rbuf_type& rb) :
            m_rbuf(&rb)
        {}
        void attach(rbuf_type& rb) { m_rbuf = &rb; }

        //--------------------------------------------------------------------
        template<class PixFmt>
        bool attach(PixFmt& pixf, int x1, int y1, int x2, int y2)
        {
            rect_i r(x1, y1, x2, y2);
            if(r.clip(rect_i(0, 0, pixf.width()-1, pixf.height()-1)))
            {
                int stride = pixf.stride();
                m_rbuf->attach(pixf.pix_ptr(r.x1, stride < 0 ? r.y2 : r.y1), 
                               (r.x2 - r.x1) + 1,
                               (r.y2 - r.y1) + 1,
                               stride);
                return true;
            }
            return false;
        }

        //--------------------------------------------------------------------
        Blender& blender() { return m_blender; }

        //--------------------------------------------------------------------
        AGG_INLINE unsigned width()  const { return m_rbuf->width();  }
        AGG_INLINE unsigned height() const { return m_rbuf->height(); }
        AGG_INLINE int      stride() const { return m_rbuf->stride(); }

        //--------------------------------------------------------------------
        AGG_INLINE       int8u* row_ptr(int y)       { return m_rbuf->row_ptr(y); }
        AGG_INLINE const int8u* row_ptr(int y) const { return m_rbuf->row_ptr(y); }
        AGG_INLINE row_data     row(int y)     const { return m_rbuf->row(y); }

        //--------------------------------------------------------------------
        AGG_INLINE int8u* pix_ptr(int x, int y) 
        { 
            return m_rbuf->row_ptr(y) + x * pix_width; 
        }

        AGG_INLINE const int8u* pix_ptr(int x, int y) const 
        { 
            return m_rbuf->row_ptr(y) + x * pix_width; 
        }

        //--------------------------------------------------------------------
        AGG_INLINE static void make_pix(int8u* p, const color_type& c)
        {
            ((value_type*)p)[order_type::R] = c.r;
            ((value_type*)p)[order_type::G] = c.g;
            ((value_type*)p)[order_type::B] = c.b;
        }

        //--------------------------------------------------------------------
        AGG_INLINE color_type pixel(int x, int y) const
        {
            value_type* p = (value_type*)m_rbuf->row_ptr(y) + x * Step;
            return color_type(p[order_type::R], 
                              p[order_type::G], 
                              p[order_type::B]);
        }

        //--------------------------------------------------------------------
        AGG_INLINE void copy_pixel(int x, int y, const color_type& c)
        {
            value_type* p = (value_type*)m_rbuf->row_ptr(x, y, 1) + x * Step;
            p[order_type::R] = c.r;
            p[order_type::G] = c.g;
            p[order_type::B] = c.b;
        }

        //--------------------------------------------------------------------
        AGG_INLINE void blend_pixel(int x, int y, const color_type& c, int8u cover)
        {
            copy_or_blend_pix((value_type*)m_rbuf->row_ptr(x, y, 1) + x * Step, c, cover);
        }


        //--------------------------------------------------------------------
        AGG_INLINE void copy_hline(int x, int y, 
                                   unsigned len, 
                                   const color_type& c)
        {
            value_type* p = (value_type*)m_rbuf->row_ptr(x, y, len) + x * Step;
            do
            {
                p[order_type::R] = c.r; 
                p[order_type::G] = c.g; 
                p[order_type::B] = c.b;
                p += Step;
            }
            while(--len);
        }


        //--------------------------------------------------------------------
        AGG_INLINE void copy_vline(int x, int y,
                                   unsigned len, 
                                   const color_type& c)
        {
            do
            {
                value_type* p = (value_type*)
                    m_rbuf->row_ptr(x, y++, 1) + x * Step;
                p[order_type::R] = c.r; 
                p[order_type::G] = c.g; 
                p[order_type::B] = c.b;
            }
            while(--len);
        }


        //--------------------------------------------------------------------
        void blend_hline(int x, int y,
                         unsigned len, 
                         const color_type& c,
                         int8u cover)
        {
            if (c.a)
            {
                value_type* p = (value_type*)
                    m_rbuf->row_ptr(x, y, len) + x * Step;

                if(c.a == base_mask && cover == cover_mask)
                {
                    do
                    {
                        p[order_type::R] = c.r; 
                        p[order_type::G] = c.g; 
                        p[order_type::B] = c.b;
                        p += Step;
                    }
                    while(--len);
                }
                else
                {
                    do
                    {
                        m_blender.blend_pix(p, c.r, c.g, c.b, c.a, cover);
                        p += Step;
                    }
                    while(--len);
                }
            }
        }


        //--------------------------------------------------------------------
        void blend_vline(int x, int y,
                         unsigned len, 
                         const color_type& c,
                         int8u cover)
        {
            if (c.a)
            {
                value_type* p;
                if(c.a == base_mask && cover == cover_mask)
                {
                    do
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x * Step;

                        p[order_type::R] = c.r; 
                        p[order_type::G] = c.g; 
                        p[order_type::B] = c.b;
                    }
                    while(--len);
                }
                else
                {
                    do
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x * Step;

                        m_blender.blend_pix(p, c.r, c.g, c.b, c.a, cover);
                    }
                    while(--len);
                }
            }
        }


        //--------------------------------------------------------------------
        void blend_solid_hspan(int x, int y,
                               unsigned len, 
                               const color_type& c,
                               const int8u* covers)
        {
            if (c.a)
            {
                value_type* p = (value_type*)
                    m_rbuf->row_ptr(x, y, len) + x * Step;

                do 
                {
                    if(c.a == base_mask && *covers == cover_mask)
                    {
                        p[order_type::R] = c.r;
                        p[order_type::G] = c.g;
                        p[order_type::B] = c.b;
                    }
                    else
                    {
                        m_blender.blend_pix(p, c.r, c.g, c.b, c.a, *covers);
                    }
                    p += Step;
                    ++covers;
                }
                while(--len);
            }
        }


        //--------------------------------------------------------------------
        void blend_solid_vspan(int x, int y,
                               unsigned len, 
                               const color_type& c,
                               const int8u* covers)
        {
            if (c.a)
            {
                do 
                {
                    value_type* p = (value_type*)
                        m_rbuf->row_ptr(x, y++, 1) + x * Step;

                    if(c.a == base_mask && *covers == cover_mask)
                    {
                        p[order_type::R] = c.r;
                        p[order_type::G] = c.g;
                        p[order_type::B] = c.b;
                    }
                    else
                    {
                        m_blender.blend_pix(p, c.r, c.g, c.b, c.a, *covers);
                    }
                    ++covers;
                }
                while(--len);
            }
        }


        //--------------------------------------------------------------------
        void copy_color_hspan(int x, int y,
                              unsigned len, 
                              const color_type* colors)
        {
            value_type* p = (value_type*)
                m_rbuf->row_ptr(x, y, len) + x * Step;

            do 
            {
                p[order_type::R] = colors->r;
                p[order_type::G] = colors->g;
                p[order_type::B] = colors->b;
                ++colors;
                p += Step;
            }
            while(--len);
        }


        //--------------------------------------------------------------------
        void copy_color_vspan(int x, int y,
                              unsigned len, 
                              const color_type* colors)
        {
            do 
            {
                value_type* p = (value_type*)
                    m_rbuf->row_ptr(x, y++, 1) + x * Step;
                p[order_type::R] = colors->r;
                p[order_type::G] = colors->g;
                p[order_type::B] = colors->b;
                ++colors;
            }
            while(--len);
        }


        //--------------------------------------------------------------------
        void blend_color_hspan(int x, int y,
                               unsigned len, 
                               const color_type* colors,
                               const int8u* covers,
                               int8u cover)
        {
            value_type* p = (value_type*)
                m_rbuf->row_ptr(x, y, len) + x * Step;

            if(covers)
            {
                do 
                {
                    copy_or_blend_pix(p, *colors++, *covers++);
                    p += Step;
                }
                while(--len);
            }
            else
            {
                if(cover == cover_mask)
                {
                    do 
                    {
                        copy_or_blend_pix(p, *colors++);
                        p += Step;
                    }
                    while(--len);
                }
                else
                {
                    do 
                    {
                        copy_or_blend_pix(p, *colors++, cover);
                        p += Step;
                    }
                    while(--len);
                }
            }
        }



        //--------------------------------------------------------------------
        void blend_color_vspan(int x, int y,
                               unsigned len, 
                               const color_type* colors,
                               const int8u* covers,
                               int8u cover)
        {
            value_type* p;
            if(covers)
            {
                do 
                {
                    p = (value_type*)
                        m_rbuf->row_ptr(x, y++, 1) + x * Step;

                    copy_or_blend_pix(p, *colors++, *covers++);
                }
                while(--len);
            }
            else
            {
                if(cover == cover_mask)
                {
                    do 
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x * Step;

                        copy_or_blend_pix(p, *colors++);
                    }
                    while(--len);
                }
                else
                {
                    do 
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x * Step;

                        copy_or_blend_pix(p, *colors++, cover);
                    }
                    while(--len);
                }
            }
        }

        //--------------------------------------------------------------------
        template<class Function> void for_each_pixel(Function f)
        {
            unsigned y;
            for(y = 0; y < height(); ++y)
            {
                row_data r = m_rbuf->row(y);
                if(r.ptr)
                {
                    unsigned len = r.x2 - r.x1 + 1;
                    value_type* p = (value_type*)
                        m_rbuf->row_ptr(r.x1, y, len) + r.x1 * Step;
                    do
                    {
                        f(p);
                        p += Step;
                    }
                    while(--len);
                }
            }
        }

        //--------------------------------------------------------------------
        template<class GammaLut> void apply_gamma_dir(const GammaLut& g)
        {
            for_each_pixel(apply_gamma_dir_rgb<color_type, order_type, GammaLut>(g));
        }

        //--------------------------------------------------------------------
        template<class GammaLut> void apply_gamma_inv(const GammaLut& g)
        {
            for_each_pixel(apply_gamma_inv_rgb<color_type, order_type, GammaLut>(g));
        }

        //--------------------------------------------------------------------
        template<class RenBuf2>
        void copy_from(const RenBuf2& from, 
                       int xdst, int ydst,
                       int xsrc, int ysrc,
                       unsigned len)
        {
            const int8u* p = from.row_ptr(ysrc);
            if(p)
            {
                memmove(m_rbuf->row_ptr(xdst, ydst, len) + xdst * pix_width, 
                        p + xsrc * pix_width, 
                        len * pix_width);
            }
        }


        //--------------------------------------------------------------------
        template<class SrcPixelFormatRenderer>
        void blend_from(const SrcPixelFormatRenderer& from, 
                        int xdst, int ydst,
                        int xsrc, int ysrc,
                        unsigned len,
                        int8u cover)
        {
            typedef typename SrcPixelFormatRenderer::order_type src_order;

            const value_type* psrc = (const value_type*)from.row_ptr(ysrc);
            if(psrc)
            {
                psrc += xsrc * 4;
                value_type* pdst = 
                    (value_type*)m_rbuf->row_ptr(xdst, ydst, len) + xdst * Step;   

                if(cover == cover_mask)
                {
                    do 
                    {
                        value_type alpha = psrc[src_order::A];
                        if(alpha)
                        {
                            if(alpha == base_mask)
                            {
                                pdst[order_type::R] = psrc[src_order::R];
                                pdst[order_type::G] = psrc[src_order::G];
                                pdst[order_type::B] = psrc[src_order::B];
                            }
                            else
                            {
                                m_blender.blend_pix(pdst, 
                                                    psrc[src_order::R],
                                                    psrc[src_order::G],
                                                    psrc[src_order::B],
                                                    alpha);
                            }
                        }
                        psrc += 4;
                        pdst += Step;
                    }
                    while(--len);
                }
                else
                {
                    color_type color;
                    do 
                    {
                        color.r = psrc[src_order::R];
                        color.g = psrc[src_order::G];
                        color.b = psrc[src_order::B];
                        color.a = psrc[src_order::A];
                        copy_or_blend_pix(pdst, color, cover);
                        psrc += 4;
                        pdst += Step;
                    }
                    while(--len);
                }
            }
        }

        //--------------------------------------------------------------------
        template<class SrcPixelFormatRenderer>
        void blend_from_color(const SrcPixelFormatRenderer& from, 
                              const color_type& color,
                              int xdst, int ydst,
                              int xsrc, int ysrc,
                              unsigned len,
                              int8u cover)
        {
            typedef typename SrcPixelFormatRenderer::value_type src_value_type;
            typedef typename SrcPixelFormatRenderer::color_type src_color_type;
            const src_value_type* psrc = (src_value_type*)from.row_ptr(ysrc);
            if(psrc)
            {
                psrc += xsrc * SrcPixelFormatRenderer::pix_step + SrcPixelFormatRenderer::pix_offset;
                value_type* pdst = 
                    (value_type*)m_rbuf->row_ptr(xdst, ydst, len) + xdst * Step;
                do 
                {
                    copy_or_blend_pix(pdst, 
                                      color, 
                                      src_color_type::int_mult_cover(*psrc, cover) >>
                                      (src_color_type::base_shift - 8));
                    psrc += SrcPixelFormatRenderer::pix_step;
                    pdst += Step;
                }
                while(--len);
            }
        }

        //--------------------------------------------------------------------
        template<class SrcPixelFormatRenderer>
        void blend_from_lut(const SrcPixelFormatRenderer& from, 
                            const color_type* color_lut,
                            int xdst, int ydst,
                            int xsrc, int ysrc,
                            unsigned len,
                            int8u cover)
        {
            typedef typename SrcPixelFormatRenderer::value_type src_value_type;
            const src_value_type* psrc = (src_value_type*)from.row_ptr(ysrc);
            if(psrc)
            {
                psrc += xsrc * SrcPixelFormatRenderer::pix_step + SrcPixelFormatRenderer::pix_offset;
                value_type* pdst = 
                    (value_type*)m_rbuf->row_ptr(xdst, ydst, len) + xdst * Step;

                if(cover == cover_mask)
                {
                    do 
                    {
                        const color_type& color = color_lut[*psrc];
                        m_blender.blend_pix(pdst, 
                                            color.r, color.g, color.b, color.a);
                        psrc += SrcPixelFormatRenderer::pix_step;
                        pdst += Step;
                    }
                    while(--len);
                }
                else
                {
                    do 
                    {
                        copy_or_blend_pix(pdst, color_lut[*psrc], cover);
                        psrc += SrcPixelFormatRenderer::pix_step;
                        pdst += Step;
                    }
                    while(--len);
                }
            }
        }

    private:
        rbuf_type* m_rbuf;
        Blender    m_blender;
    };

    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba8,  order_rgb>, rendering_buffer> pixfmt_rgb24;    //----pixfmt_rgb24
    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba8,  order_bgr>, rendering_buffer> pixfmt_bgr24;    //----pixfmt_bgr24
    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba16, order_rgb>, rendering_buffer> pixfmt_rgb48;    //----pixfmt_rgb48
    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba16, order_bgr>, rendering_buffer> pixfmt_bgr48;    //----pixfmt_bgr48

    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba8,  order_rgb>, rendering_buffer> pixfmt_rgb24_pre; //----pixfmt_rgb24_pre
    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba8,  order_bgr>, rendering_buffer> pixfmt_bgr24_pre; //----pixfmt_bgr24_pre
    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba16, order_rgb>, rendering_buffer> pixfmt_rgb48_pre; //----pixfmt_rgb48_pre
    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba16, order_bgr>, rendering_buffer> pixfmt_bgr48_pre; //----pixfmt_bgr48_pre

    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba8,  order_rgb>, rendering_buffer, 4> pixfmt_rgbx32;    //----pixfmt_rgb24
    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba8,  order_bgr>, rendering_buffer, 4> pixfmt_bgrx32;    //----pixfmt_bgr24
    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba16, order_rgb>, rendering_buffer, 4> pixfmt_rgbx64;    //----pixfmt_rgb48
    typedef pixfmt_alpha_blend_rgb<blender_rgb<rgba16, order_bgr>, rendering_buffer, 4> pixfmt_bgrx64;    //----pixfmt_bgr48
    
    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba8,  order_rgb>, rendering_buffer, 4> pixfmt_rgbx32_pre; //----pixfmt_rgb24_pre
    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba8,  order_bgr>, rendering_buffer, 4> pixfmt_bgrx32_pre; //----pixfmt_bgr24_pre
    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba16, order_rgb>, rendering_buffer, 4> pixfmt_rgbx64_pre; //----pixfmt_rgb48_pre
    typedef pixfmt_alpha_blend_rgb<blender_rgb_pre<rgba16, order_bgr>, rendering_buffer, 4> pixfmt_bgrx64_pre; //----pixfmt_bgr48_pre
    
    //-----------------------------------------------------pixfmt_rgb24_gamma
    template<class Gamma> class pixfmt_rgb24_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_rgb, Gamma>, rendering_buffer>
    {
    public:
        pixfmt_rgb24_gamma(rendering_buffer& rb, const Gamma& g) :
            pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_rgb, Gamma>, rendering_buffer>(rb) 
        {
            this->blender().gamma(g);
        }
    };
        
    //-----------------------------------------------------pixfmt_bgr24_gamma
    template<class Gamma> class pixfmt_bgr24_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_bgr, Gamma>, rendering_buffer>
    {
    public:
        pixfmt_bgr24_gamma(rendering_buffer& rb, const Gamma& g) :
            pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_bgr, Gamma>, rendering_buffer>(rb) 
        {
            this->blender().gamma(g);
        }
    };

    //-----------------------------------------------------pixfmt_rgb48_gamma
    template<class Gamma> class pixfmt_rgb48_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_rgb, Gamma>, rendering_buffer>
    {
    public:
        pixfmt_rgb48_gamma(rendering_buffer& rb, const Gamma& g) :
            pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_rgb, Gamma>, rendering_buffer>(rb) 
        {
            this->blender().gamma(g);
        }
    };
        
    //-----------------------------------------------------pixfmt_bgr48_gamma
    template<class Gamma> class pixfmt_bgr48_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_bgr, Gamma>, rendering_buffer>
    {
    public:
        pixfmt_bgr48_gamma(rendering_buffer& rb, const Gamma& g) :
            pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_bgr, Gamma>, rendering_buffer>(rb) 
        {
            this->blender().gamma(g);
        }
    };

    //-----------------------------------------------------pixfmt_rgbx32_gamma
    template<class Gamma> class pixfmt_rgbx32_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_rgb, Gamma>, rendering_buffer, 4>
    {
    public:
        pixfmt_rgbx32_gamma(rendering_buffer& rb, const Gamma& g) :
        pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_rgb, Gamma>, rendering_buffer, 4>(rb) 
        {
            this->blender().gamma(g);
        }
    };
    
    //-----------------------------------------------------pixfmt_bgrx32_gamma
    template<class Gamma> class pixfmt_bgrx32_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_bgr, Gamma>, rendering_buffer, 4>
    {
    public:
        pixfmt_bgrx32_gamma(rendering_buffer& rb, const Gamma& g) :
        pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba8, order_bgr, Gamma>, rendering_buffer, 4>(rb) 
        {
            this->blender().gamma(g);
        }
    };
    
    //-----------------------------------------------------pixfmt_rgbx64_gamma
    template<class Gamma> class pixfmt_rgbx64_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_rgb, Gamma>, rendering_buffer, 4>
    {
    public:
        pixfmt_rgbx64_gamma(rendering_buffer& rb, const Gamma& g) :
        pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_rgb, Gamma>, rendering_buffer, 4>(rb) 
        {
            this->blender().gamma(g);
        }
    };
    
    //-----------------------------------------------------pixfmt_bgrx64_gamma
    template<class Gamma> class pixfmt_bgrx64_gamma : 
    public pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_bgr, Gamma>, rendering_buffer, 4>
    {
    public:
        pixfmt_bgrx64_gamma(rendering_buffer& rb, const Gamma& g) :
        pixfmt_alpha_blend_rgb<blender_rgb_gamma<rgba16, order_bgr, Gamma>, rendering_buffer, 4>(rb) 
        {
            this->blender().gamma(g);
        }
    };
    
    
}

#endif

