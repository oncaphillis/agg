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

#ifndef AGG_PIXFMT_DEVICE_N_INCLUDED
#define AGG_PIXFMT_DEVICE_N_INCLUDED

#include <string.h>
#include "agg_basics.h"
#include "agg_color_device_na.h"
#include "agg_rendering_buffer.h"

namespace agg
{

    //=====================================================apply_gamma_dir_device_n
    template<class ColorT, class Order, class GammaLut>
    class apply_gamma_dir_device_n
    {
    public:
        typedef typename ColorT::value_type value_type;

        apply_gamma_dir_device_n(const GammaLut& gamma) : m_gamma(gamma) {}

        AGG_INLINE void operator () (value_type * p)
        {
            p[Order::C] = m_gamma.dir(p[Order::C]);
            p[Order::M] = m_gamma.dir(p[Order::M]);
            p[Order::Y] = m_gamma.dir(p[Order::Y]);
            p[Order::K] = m_gamma.dir(p[Order::K]);
        }

    private:
        const GammaLut& m_gamma;
    };



    //=====================================================apply_gamma_inv_cmyk
    template<class ColorT, class Order, class GammaLut>
    class apply_gamma_inv_device_n
    {
    public:
        typedef typename ColorT::value_type value_type;

        apply_gamma_inv_device_n(const GammaLut& gamma) : m_gamma(gamma) {}

        AGG_INLINE void operator () (value_type* p)
        {
            p[Order::C] = m_gamma.inv(p[Order::C]);
            p[Order::M] = m_gamma.inv(p[Order::M]);
            p[Order::Y] = m_gamma.inv(p[Order::Y]);
            p[Order::K] = m_gamma.inv(p[Order::K]);
        }

    private:
        const GammaLut& m_gamma;
    };


    //=========================================================blender_device_n
    template<class ColorT, class Order>
    struct blender_device_n
    {
        typedef ColorT color_type;
        typedef Order order_type;
        typedef typename color_type::value_type value_type;
        typedef typename color_type::calc_type calc_type;
        enum base_scale_e { base_shift = color_type::base_shift };

        //--------------------------------------------------------------------
        static AGG_INLINE void blend_pix(value_type* p, 
                                         unsigned cc, unsigned cm, unsigned cy, unsigned ck, 
                                         unsigned alpha, 
                                         unsigned cover=0)
        {
            p[Order::C] += (value_type)(((cc - p[Order::C]) * alpha) >> base_shift);
            p[Order::M] += (value_type)(((cm - p[Order::M]) * alpha) >> base_shift);
            p[Order::Y] += (value_type)(((cy - p[Order::Y]) * alpha) >> base_shift);
            p[Order::K] += (value_type)(((ck - p[Order::K]) * alpha) >> base_shift);
        }
    };


    //======================================================blender_device_n_pre
    template<class ColorT, class Order>
    struct blender_device_n_pre
    {
        typedef ColorT color_type;
        typedef Order order_type;
        typedef typename color_type::value_type value_type;
        typedef typename color_type::calc_type calc_type;
        enum base_scale_e { base_shift = color_type::base_shift };

        //--------------------------------------------------------------------
        static AGG_INLINE void blend_pix(value_type* p, 
                                         unsigned cc, unsigned cm, unsigned cy, unsigned ck,
                                         unsigned alpha,
                                         unsigned cover)
        {
            alpha = color_type::base_mask - alpha;
            cover = (cover + 1) << (base_shift - 8);

            p[Order::C] = (value_type)((p[Order::C] * alpha + cc * cover) >> base_shift);
            p[Order::M] = (value_type)((p[Order::M] * alpha + cm * cover) >> base_shift);
            p[Order::Y] = (value_type)((p[Order::Y] * alpha + cy * cover) >> base_shift);
            p[Order::K] = (value_type)((p[Order::K] * alpha + ck * cover) >> base_shift);
        }

        //--------------------------------------------------------------------
        static AGG_INLINE void blend_pix(value_type* p, 
                                         unsigned cc, unsigned cm, unsigned cy, unsigned ck,
                                         unsigned alpha)
        {
            alpha = color_type::base_mask - alpha;
            p[Order::C] = (value_type)(((p[Order::C] * alpha) >> base_shift) + cc);
            p[Order::M] = (value_type)(((p[Order::M] * alpha) >> base_shift) + cm);
            p[Order::Y] = (value_type)(((p[Order::Y] * alpha) >> base_shift) + cy);
            p[Order::K] = (value_type)(((p[Order::K] * alpha) >> base_shift) + ck);
        }
    };



    //===================================================blender_device_n_gamma
    template<class ColorT, class Order, class Gamma,int N>
    class blender_device_n_gamma
    {
    public:
        typedef ColorT color_type;
        typedef Order order_type;
        typedef Gamma gamma_type;
        typedef typename color_type::value_type value_type;
        typedef typename color_type::calc_type calc_type;
        enum base_scale_e { base_shift = color_type::base_shift };

        //--------------------------------------------------------------------
        blender_device_n_gamma() : m_gamma(0)
        {
        }

        void gamma(const gamma_type& g) { m_gamma = &g; }

        //--------------------------------------------------------------------
        AGG_INLINE void blend_pix(value_type* p, 
                                  unsigned cc, unsigned cm, unsigned cy, unsigned ck,
                                  unsigned alpha, 
                                  unsigned cover=0)
        {
            calc_type c = m_gamma->dir(p[Order::C]);
            calc_type m = m_gamma->dir(p[Order::M]);
            calc_type y = m_gamma->dir(p[Order::Y]);
            calc_type k = m_gamma->dir(p[Order::K]);

            p[Order::C] = m_gamma->inv((((m_gamma->dir(cc) - c) * alpha) >> base_shift) + c);
            p[Order::M] = m_gamma->inv((((m_gamma->dir(cm) - m) * alpha) >> base_shift) + m);
            p[Order::Y] = m_gamma->inv((((m_gamma->dir(cy) - y) * alpha) >> base_shift) + y);
            p[Order::K] = m_gamma->inv((((m_gamma->dir(ck) - k) * alpha) >> base_shift) + k);
        }

    private:
        const gamma_type* m_gamma;
    };



    
    //==================================================pixfmt_alpha_device_n_rgb
    template<class Blender, class RenBuf>
    class pixfmt_alpha_blend_device_n
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
            pix_width  = sizeof(value_type) * 4
        };

    private:
        //--------------------------------------------------------------------
        AGG_INLINE void copy_or_blend_pix(value_type* p, 
                                          const color_type& c, 
                                          unsigned cover)
        {
            if (c.a)
            {
                calc_type alpha = (calc_type(c.a) * (cover + 1)) >> 8;
                if(alpha == base_mask)
                {
                    p[order_type::C] = c.c;
                    p[order_type::M] = c.m;
                    p[order_type::Y] = c.y;
                    p[order_type::K] = c.k;
                }
                else
                {
                    m_blender.blend_pix(p, c.c, c.m, c.y, c.k, alpha, cover);
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
                    p[order_type::C] = c.c;
                    p[order_type::M] = c.m;
                    p[order_type::Y] = c.y;
                    p[order_type::K] = c.k;
                }
                else
                {
                    m_blender.blend_pix(p, c.c, c.m, c.y, c.k, c.a);
                }
            }
        }


    public:
        //--------------------------------------------------------------------
        explicit pixfmt_alpha_blend_device_n(rbuf_type& rb) :
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
        AGG_INLINE static void make_pix(int8u* p, const color_type& cl)
        {
            ((value_type*)p)[order_type::C] = cl.c;
            ((value_type*)p)[order_type::M] = cl.m;
            ((value_type*)p)[order_type::Y] = cl.y;
            ((value_type*)p)[order_type::K] = cl.k;
        }

        //--------------------------------------------------------------------
        AGG_INLINE color_type pixel(int x, int y) const
        {
            value_type* p = (value_type*)m_rbuf->row_ptr(y) + x + x + x + x;
            return color_type(p[order_type::C], 
                              p[order_type::M], 
                              p[order_type::Y],
                              p[order_type::K]);
        }

        //--------------------------------------------------------------------
        AGG_INLINE void copy_pixel(int x, int y, const color_type& cl)
        {
            value_type* p = (value_type*)m_rbuf->row_ptr(x, y, 1) + x + x + x + x;
            p[order_type::C] = cl.c;
            p[order_type::M] = cl.m;
            p[order_type::Y] = cl.y;
            p[order_type::K] = cl.k;
        }

        //--------------------------------------------------------------------
        AGG_INLINE void blend_pixel(int x, int y, const color_type& c, int8u cover)
        {
            copy_or_blend_pix((value_type*)m_rbuf->row_ptr(x, y, 1) + x + x + x +x, c, cover);
        }


        //--------------------------------------------------------------------
        AGG_INLINE void copy_hline(int x, int y, 
                                   unsigned len, 
                                   const color_type& cl)
        {
            value_type* p = (value_type*)m_rbuf->row_ptr(x, y, len) + x + x + x + x;
            do
            {
                p[order_type::C] = cl.c; 
                p[order_type::M] = cl.m; 
                p[order_type::Y] = cl.y;
                p[order_type::K] = cl.k;
                p += 4;
            }
            while(--len);
        }


        //--------------------------------------------------------------------
        AGG_INLINE void copy_vline(int x, int y,
                                   unsigned len, 
                                   const color_type& cl)
        {
            do
            {
                value_type* p = (value_type*)
                    m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;
                p[order_type::C] = cl.c;
                p[order_type::M] = cl.m; 
                p[order_type::Y] = cl.y;
                p[order_type::K] = cl.k;
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
                    m_rbuf->row_ptr(x, y, len) + x + x + x + x;

                calc_type alpha = (calc_type(c.a) * (calc_type(cover) + 1)) >> 8;
                if(alpha == base_mask)
                {
                    do
                    {
                        p[order_type::C] = c.c; 
                        p[order_type::M] = c.m; 
                        p[order_type::Y] = c.y;
                        p[order_type::K] = c.k;
                        p += 4;
                    }
                    while(--len);
                }
                else
                {
                    do
                    {
                        m_blender.blend_pix(p, c.c, c.m, c.y, c.k, alpha, cover);
                        p += 4;
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
                calc_type alpha = (calc_type(c.a) * (cover + 1)) >> 8;
                if(alpha == base_mask)
                {
                    do
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;

                        p[order_type::C] = c.c; 
                        p[order_type::M] = c.m; 
                        p[order_type::Y] = c.y;
                        p[order_type::K] = c.k;
                    }
                    while(--len);
                }
                else
                {
                    do
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;

                        m_blender.blend_pix(p, c.c, c.m, c.y, c.k, alpha, cover);
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
                    m_rbuf->row_ptr(x, y, len) + x + x + x + x;

                do 
                {
                    calc_type alpha = (calc_type(c.a) * (calc_type(*covers) + 1)) >> 8;
                    if(alpha == base_mask)
                    {
                        p[order_type::C] = c.c;
                        p[order_type::M] = c.m;
                        p[order_type::Y] = c.y;
                        p[order_type::K] = c.k;
                    }
                    else
                    {
                        m_blender.blend_pix(p, c.c, c.m, c.y, c.k, alpha, *covers);
                    }
                    p += 4;
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
                        m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;

                    calc_type alpha = (calc_type(c.a) * (calc_type(*covers) + 1)) >> 8;
                    if(alpha == base_mask)
                    {
                        p[order_type::C] = c.c;
                        p[order_type::M] = c.m;
                        p[order_type::Y] = c.y;
                        p[order_type::K] = c.k;
                    }
                    else
                    {
                        m_blender.blend_pix(p, c.c, c.m, c.y, c.k, alpha, *covers);
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
                m_rbuf->row_ptr(x, y, len) + x + x + x + x;

            do 
            {
                p[order_type::C] = colors->c;
                p[order_type::M] = colors->m;
                p[order_type::Y] = colors->y;
                p[order_type::K] = colors->k;
                ++colors;
                p += 4;
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
                    m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;
                p[order_type::C] = colors->c;
                p[order_type::M] = colors->m;
                p[order_type::Y] = colors->y;
                p[order_type::K] = colors->k;
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
                m_rbuf->row_ptr(x, y, len) + x + x + x + x;

            if(covers)
            {
                do 
                {
                    copy_or_blend_pix(p, *colors++, *covers++);
                    p += 4;
                }
                while(--len);
            }
            else
            {
                if(cover == 255)
                {
                    do 
                    {
                        copy_or_blend_pix(p, *colors++);
                        p += 4;
                    }
                    while(--len);
                }
                else
                {
                    do 
                    {
                        copy_or_blend_pix(p, *colors++, cover);
                        p += 4;
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
                        m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;

                    copy_or_blend_pix(p, *colors++, *covers++);
                }
                while(--len);
            }
            else
            {
                if(cover == 255)
                {
                    do 
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;

                        copy_or_blend_pix(p, *colors++);
                    }
                    while(--len);
                }
                else
                {
                    do 
                    {
                        p = (value_type*)
                            m_rbuf->row_ptr(x, y++, 1) + x + x + x + x;

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
                        m_rbuf->row_ptr(r.x1, y, len) + r.x1 * 4;
                    do
                    {
                        f(p);
                        p += 4;
                    }
                    while(--len);
                }
            }
        }

        //--------------------------------------------------------------------
        template<class GammaLut> void apply_gamma_dir(const GammaLut& g)
        {
            for_each_pixel(apply_gamma_dir_cmyk<color_type, order_type, GammaLut>(g));
        }

        //--------------------------------------------------------------------
        template<class GammaLut> void apply_gamma_inv(const GammaLut& g)
        {
            for_each_pixel(apply_gamma_inv_cmyk<color_type, order_type, GammaLut>(g));
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
                psrc += xsrc * 5;
                value_type* pdst = 
                    (value_type*)m_rbuf->row_ptr(xdst, ydst, len) + xdst * 4;   

                if(cover == 255)
                {
                    do 
                    {
                        value_type alpha = psrc[src_order::A];
                        if(alpha)
                        {
                            if(alpha == base_mask)
                            {
                                pdst[order_type::C] = psrc[src_order::C];
                                pdst[order_type::M] = psrc[src_order::M];
                                pdst[order_type::Y] = psrc[src_order::Y];
                                pdst[order_type::K] = psrc[src_order::K];
                            }
                            else
                            {
                                m_blender.blend_pix(pdst, 
                                                    psrc[src_order::C],
                                                    psrc[src_order::M],
                                                    psrc[src_order::Y],
                                                    psrc[src_order::K],
                                                    alpha);
                            }
                        }
                        psrc += 4;
                        pdst += 5;
                    }
                    while(--len);
                }
                else
                {
                    color_type color;
                    do 
                    {
                        color.c = psrc[src_order::C];
                        color.m = psrc[src_order::M];
                        color.y = psrc[src_order::Y];
                        color.k = psrc[src_order::K];
                        color.a = psrc[src_order::A];
                        copy_or_blend_pix(pdst, color, cover);
                        psrc += 5;
                        pdst += 4;
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
            const src_value_type* psrc = (src_value_type*)from.row_ptr(ysrc);
            if(psrc)
            {
                value_type* pdst = 
                    (value_type*)m_rbuf->row_ptr(xdst, ydst, len) + xdst * 4;
                do 
                {
                    copy_or_blend_pix(pdst, 
                                      color, 
                                      (*psrc * cover + base_mask) >> base_shift);
                    ++psrc;
                    pdst += 4;
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
                value_type* pdst = 
                    (value_type*)m_rbuf->row_ptr(xdst, ydst, len) + xdst * 4;

                if(cover == 255)
                {
                    do 
                    {
                        const color_type& color = color_lut[*psrc];
                        m_blender.blend_pix(pdst, 
                                            color.r, color.g, color.b, color.a);
                        ++psrc;
                        pdst += 4;
                    }
                    while(--len);
                }
                else
                {
                    do 
                    {
                        copy_or_blend_pix(pdst, color_lut[*psrc], cover);
                        ++psrc;
                        pdst += 4;
                    }
                    while(--len);
                }
            }
        }

    private:
        rbuf_type* m_rbuf;
        Blender    m_blender;
    };

    //----pixfmt_device_n_32
    typedef pixfmt_alpha_blend_device_n<blender_device_n<device_na8<4>,  order_cmyk>, rendering_buffer>   pixfmt_device_n_32;
}

#endif // AGG_PIXFMT_DEVICE_N_INCLUDED

