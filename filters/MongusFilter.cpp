/******************************************************************************
 * Copyright (c) 2016-2017, Bradley J Chambers (brad.chambers@gmail.com)
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 ****************************************************************************/

#include "MongusFilter.hpp"

#include <io/BufferReader.hpp>
#include <pdal/EigenUtils.hpp>
#include <pdal/PipelineManager.hpp>
#include <pdal/Segmentation.hpp>
#include <pdal/pdal_macros.hpp>
#include <pdal/util/FileUtils.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include "private/DimRange.hpp"

namespace pdal
{

static PluginInfo const s_info =
    PluginInfo("filters.mongus", "Mongus and Zalik (2012)",
               "http://pdal.io/stages/filters.mongus.html");

CREATE_STATIC_PLUGIN(1, 0, MongusFilter, Filter, s_info)

std::string MongusFilter::getName() const
{
    return s_info.name;
}

void MongusFilter::addArgs(ProgramArgs& args)
{
    args.add("cell", "Cell size", m_cellSize, 1.0);
    args.add("k", "Stdev multiplier for threshold", m_k, 3.0);
    args.add("l", "Max level", m_l, 8);
    args.add("ignore", "Ignore values", m_ignored);
    args.add("last", "Consider last returns only?", m_lastOnly, true);
    args.add("dir", "Optional output directory for debugging", m_dir);
}

void MongusFilter::addDimensions(PointLayoutPtr layout)
{
    layout->registerDim(Dimension::Id::Classification);
}

void MongusFilter::ready(PointTableRef table)
{
    if (m_dir.empty())
        return;

    if (!FileUtils::directoryExists(m_dir))
        throwError("Output directory '" + m_dir + "' does not exist");
}

size_t MongusFilter::getColIndex(double x, double cell_size)
{
    return static_cast<size_t>(floor((x - m_bounds.minx) / cell_size));
}

size_t MongusFilter::getRowIndex(double y, double cell_size)
{
    return static_cast<size_t>(floor((y - m_bounds.miny) / cell_size));
}

void MongusFilter::writeControl(std::vector<double> cx, std::vector<double> cy,
                                std::vector<double> cz, std::string filename)
{
    if (m_dir.empty())
        return;

    PipelineManager m;

    PointTable table;
    PointViewPtr view(new PointView(table));

    table.layout()->registerDim(Dimension::Id::X);
    table.layout()->registerDim(Dimension::Id::Y);
    table.layout()->registerDim(Dimension::Id::Z);

    PointId i = 0;
    for (size_t j = 0; j < cz.size(); ++j)
    {
        if (std::isnan(cx[j]) || std::isnan(cy[j]) || std::isnan(cz[j]))
            continue;
        view->setField(Dimension::Id::X, i, cx[j]);
        view->setField(Dimension::Id::Y, i, cy[j]);
        view->setField(Dimension::Id::Z, i, cz[j]);
        i++;
    }

    BufferReader r;
    r.addView(view);

    std::string fname = FileUtils::toAbsolutePath(filename, m_dir);
    Stage& w = m.makeWriter(fname, "writers.las", r);
    w.prepare(table);
    w.execute(table);
}

std::vector<PointId> MongusFilter::processGround(PointViewPtr view)
{
    std::vector<PointId> groundIdx;

    view->calculateBounds(m_bounds);

    m_cols = ((m_bounds.maxx - m_bounds.minx) / m_cellSize) + 1;
    m_rows = ((m_bounds.maxy - m_bounds.miny) / m_cellSize) + 1;

    // create control points matrix at default cell size
    std::vector<double> cx(m_rows * m_cols,
                           std::numeric_limits<double>::quiet_NaN());

    std::vector<double> cy(m_rows * m_cols,
                           std::numeric_limits<double>::quiet_NaN());

    std::vector<double> cz(m_rows * m_cols,
                           std::numeric_limits<double>::quiet_NaN());

    // find initial set of Z minimums at native resolution
    for (PointId i = 0; i < view->size(); ++i)
    {
        double x = view->getFieldAs<double>(Dimension::Id::X, i);
        double y = view->getFieldAs<double>(Dimension::Id::Y, i);
        double z = view->getFieldAs<double>(Dimension::Id::Z, i);

        size_t c = getColIndex(x, m_cellSize);
        size_t r = getRowIndex(y, m_cellSize);

        if (z < cz[c * m_rows + r] || std::isnan(cz[c * m_rows + r]))
        {
            cx[c * m_rows + r] = x;
            cy[c * m_rows + r] = y;
            cz[c * m_rows + r] = z;
        }
    }

    // In our case, 2D structural elements of circular shape are employed and
    // sufficient accuracy is achieved by using a larger window size for opening
    // (W11) than for closing (W9).
    std::vector<double> erode = eigen::erodeDiamond(cz, m_rows, m_cols, 11);
    std::vector<double> mo = eigen::dilateDiamond(erode, m_rows, m_cols, 11);
    writeControl(cx, cy, mo, "grid_open.laz");
    std::vector<double> dilate = eigen::dilateDiamond(mo, m_rows, m_cols, 9);
    std::vector<double> mc = eigen::erodeDiamond(dilate, m_rows, m_cols, 9);
    writeControl(cx, cy, mc, "grid_close.laz");

    // ...in order to minimize the distortions caused by such filtering, the
    // output points ... are compared to C and only ci with significantly lower
    // elevation [are] replaced... In our case, d = 1.0 m was used.
    for (size_t i = 0; i < cz.size(); ++i)
    {
        if ((mc[i] - cz[i]) >= 1.0)
            cz[i] = mc[i];
    }
    // cz is still at native resolution, with low points replaced by
    // morphological operators
    writeControl(cx, cy, cz, "grid_mins_adjusted.laz");

    // downsample control at max_level
    int level = m_l;
    double cur_cell_size = m_cellSize * std::pow(2, level);
    // for max level = 8 and cell size 1, this is 256

    std::vector<double> x_prev, y_prev, z_prev;

    // Top-level control samples are assumed to be ground points, no filtering
    // is applied.
    size_t nr = std::ceil(m_rows / cur_cell_size);
    size_t nc = std::ceil(m_cols / cur_cell_size);
    size_t pr = nr;
    size_t pc = nc;
    downsampleMin(&cx, &cy, &cz, &x_prev, &y_prev, &z_prev, cur_cell_size, nr,
                  nc);
    // x|y|z_prev are control points downsampled to coarsest resolution for
    // the hierarchy, e.g., for 512x512, this would be 2x2
    writeControl(x_prev, y_prev, z_prev, "control_init.laz");

    // Point-filtering is performed iteratively at each level of the
    // control-points hierarchy in a top-down fashion
    for (auto l = level - 1; l > 0; --l)
    {
        cur_cell_size /= 2;
        // 128, 64, 32, 16, 8, 4, 1

        // compute TPS with update control at level

        // The interpolated surface is estimated based on the filtered set of
        // TPS control-points at the previous level of hierarchy
        // MatrixXd surface = TPS(x_prev, y_prev, z_prev, cur_cell_size);
        // 4x4, 8x8, 16x16, 32x32, 64x64, 128x128, 256x256

        // downsample control at level
        std::vector<double> x_samp, y_samp, z_samp;
        nr = std::ceil(m_rows / cur_cell_size);
        nc = std::ceil(m_cols / cur_cell_size);
        downsampleMin(&cx, &cy, &cz, &x_samp, &y_samp, &z_samp, cur_cell_size,
                      nr, nc);
        // 4x4, 8x8, 16x16, 32x32, 64x64, 128x128, 256x256

        // If the number of indices exceeds some threshold, then use the sampled
        // spline method, rather than computing the full spline.
        std::vector<double> surface;
        if (nr * nc > 10)
            surface = eigen::computeSampledSpline(x_prev, y_prev, z_prev, pr,
                                                  pc, x_samp, y_samp, nr, nc);
        else
            surface = eigen::computeSpline(x_prev, y_prev, z_prev, pr, pc,
                                           x_samp, y_samp, nr, nc);

        char buffer[256];
        sprintf(buffer, "interp_surface_%d.laz", l);
        std::string name(buffer);
        writeControl(x_samp, y_samp, surface, name);

        std::vector<double> R;
        std::transform(z_samp.begin(), z_samp.end(), surface.begin(),
                       std::back_inserter(R),
                       [](double l, double r) { return l - r; });

        char rbuf[256];
        sprintf(rbuf, "residual_%d.laz", l);
        std::string rbufn(rbuf);
        writeControl(x_samp, y_samp, R, rbufn);

        // compute top-hat (white) of residuals
        std::vector<double> erodeR = eigen::erodeDiamond(R, nr, nc, 2 * l);
        std::vector<double> openR = eigen::dilateDiamond(erodeR, nr, nc, 2 * l);
        std::vector<double> T;
        std::transform(R.begin(), R.end(), openR.begin(), std::back_inserter(T),
                       [](double l, double r) { return l - r; });

        char mbuf[256];
        sprintf(mbuf, "median_%d.laz", l);
        std::string mbufn(mbuf);
        writeControl(x_samp, y_samp, T, mbufn);

        // compute mean and stddev of top-hat values
        double M1, M2;
        M1 = M2 = 0.0;
        point_count_t n(0);
        for (size_t i = 0; i < T.size(); ++i)
        {
            point_count_t n1(n);
            n++;
            double delta = T[i] - M1;
            double delta_n = delta / n;
            double term1 = delta * delta_n * n1;
            M1 += delta_n;
            M2 += term1;
        }
        double mean = M1;
        double stddev = std::sqrt(M2 / (n - 1.0));

        // If the TPS control-point is recognized as a non-ground point, it is
        // replaced by the interpolated point. The time complexity of the
        // approach is reduced by filtering only the control-points in each
        // iteration.
        for (size_t i = 0; i < T.size(); ++i)
        {
            if (T[i] > (mean + 3 * stddev))
                z_samp[i] = std::numeric_limits<double>::quiet_NaN();
        }

        char buf2[256];
        sprintf(buf2, "adjusted_control_%d.laz", l);
        std::string name2(buf2);
        writeControl(x_samp, y_samp, z_samp, name2);

        char buf3[256];
        sprintf(buf3, "prev_control_%d.laz", l);
        std::string name3(buf3);
        writeControl(x_prev, y_prev, z_prev, name3);

        x_prev.swap(x_samp);
        y_prev.swap(y_samp);
        z_prev.swap(z_samp);
        pr = nr;
        pc = nc;
    }

    std::vector<double> surface = eigen::computeSampledSpline(
        x_prev, y_prev, z_prev, pr, pc, cx, cy, m_rows, m_cols);

    writeControl(cx, cy, surface, "final_surface.laz");

    // apply final filtering (top hat) using raw points against TPS
    std::vector<double> R;
    std::transform(cz.begin(), cz.end(), surface.begin(), std::back_inserter(R),
                   [](double l, double r) { return l - r; });

    writeControl(cx, cy, R, "final_residual.laz");

    // compute top-hat (white) of residuals
    std::vector<double> erodeR = eigen::erodeDiamond(R, m_rows, m_cols, 1);
    std::vector<double> openR = eigen::dilateDiamond(erodeR, m_rows, m_cols, 1);
    std::vector<double> T;
    std::transform(R.begin(), R.end(), openR.begin(), std::back_inserter(T),
                   [](double l, double r) { return l - r; });

    writeControl(cx, cy, T, "final_median.laz");

    // compute mean and stddev of top-hat values
    double M1, M2;
    M1 = M2 = 0.0;
    point_count_t n(0);
    for (size_t i = 0; i < T.size(); ++i)
    {
        point_count_t n1(n);
        n++;
        double delta = T[i] - M1;
        double delta_n = delta / n;
        double term1 = delta * delta_n * n1;
        M1 += delta_n;
        M2 += term1;
    }
    double mean = M1;
    double stddev = std::sqrt(M2 / (n - 1.0));

    // If the TPS control-point is recognized as a non-ground point, it is
    // replaced by the interpolated point. The time complexity of the
    // approach is reduced by filtering only the control-points in each
    // iteration.
    for (size_t i = 0; i < T.size(); ++i)
    {
        if (T[i] > (mean + 3 * stddev))
            cz[i] = std::numeric_limits<double>::quiet_NaN();
    }

    writeControl(cx, cy, cz, "final_control.laz");

    // ...the LiDAR points are filtered only at the bottom level.
    for (PointId i = 0; i < view->size(); ++i)
    {
        using namespace Dimension;

        double x = view->getFieldAs<double>(Id::X, i);
        double y = view->getFieldAs<double>(Id::Y, i);
        double z = view->getFieldAs<double>(Id::Z, i);

        size_t c = getColIndex(x, m_cellSize);
        size_t r = getRowIndex(y, m_cellSize);

        // double res = z - surface[c * m_rows + r];
        double res = z - cz[c * m_rows + r];
        if (res < 1.0)
            groundIdx.push_back(i);
    }

    return groundIdx;
}

void MongusFilter::downsampleMin(
    std::vector<double>* cx, std::vector<double>* cy, std::vector<double>* cz,
    std::vector<double>* dcx, std::vector<double>* dcy,
    std::vector<double>* dcz, double cell_size, size_t nr, size_t nc)
{
    dcx->assign(nr * nc, std::numeric_limits<double>::quiet_NaN());
    dcy->assign(nr * nc, std::numeric_limits<double>::quiet_NaN());
    dcz->assign(nr * nc, std::numeric_limits<double>::quiet_NaN());

    for (size_t c = 0; c < m_cols; ++c)
    {
        for (size_t r = 0; r < m_rows; ++r)
        {
            if (std::isnan((*cz)[c * m_rows + r]))
                continue;

            int rr = std::floor(r / cell_size);
            int cc = std::floor(c / cell_size);

            if ((*cz)[c * m_rows + r] < (*dcz)[cc * nr + rr] ||
                std::isnan((*dcz)[cc * nr + rr]))
            {
                (*dcx)[cc * nr + rr] = (*cx)[c * m_rows + r];
                (*dcy)[cc * nr + rr] = (*cy)[c * m_rows + r];
                (*dcz)[cc * nr + rr] = (*cz)[c * m_rows + r];
            }
        }
    }
}

PointViewSet MongusFilter::run(PointViewPtr view)
{
    PointViewSet viewSet;
    if (!view->size())
        return viewSet;

    // Segment input view into ignored/kept views.
    PointViewPtr ignoredView = view->makeNew();
    PointViewPtr keptView = view->makeNew();
    if (m_ignored.m_id == Dimension::Id::Unknown)
        keptView->append(*view);
    else
        Segmentation::ignoreDimRange(m_ignored, view, keptView, ignoredView);

    // Segment kept view into last/other-than-last return views.
    PointViewPtr lastView = keptView->makeNew();
    PointViewPtr nonlastView = keptView->makeNew();
    if (m_lastOnly)
        Segmentation::segmentLastReturns(keptView, lastView, nonlastView);
    else
        lastView->append(*keptView);

    for (PointId i = 0; i < nonlastView->size(); ++i)
        nonlastView->setField(Dimension::Id::Classification, i, 1);

    std::vector<PointId> idx = processGround(view);

    PointViewPtr outView = view->makeNew();
    if (!idx.empty())
    {

        log()->get(LogLevel::Debug2)
            << "Labeled " << idx.size() << " ground returns!\n";

        // set the classification label of ground returns as 2
        // (corresponding to ASPRS LAS specification)
        for (const auto& i : idx)
        {
            lastView->setField(Dimension::Id::Classification, i, 2);
        }

        outView->append(*ignoredView);
        outView->append(*nonlastView);
        outView->append(*lastView);
    }
    else
    {
        if (idx.empty())
            log()->get(LogLevel::Debug2) << "Filtered cloud has no "
                                            "ground returns!\n";

        // return the input buffer unchanged
        outView->append(*view);
    }
    viewSet.insert(outView);

    return viewSet;
}

} // namespace pdal
