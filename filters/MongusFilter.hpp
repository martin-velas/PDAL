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

#pragma once

#include <pdal/Filter.hpp>
#include <pdal/plugin.hpp>

#include "private/DimRange.hpp"

#include <Eigen/Dense>

#include <string>
#include <vector>

extern "C" int32_t MongusFilter_ExitFunc();
extern "C" PF_ExitFunc MongusFilter_InitPlugin();

namespace pdal
{

class PointLayout;
class PointView;

class PDAL_DLL MongusFilter : public Filter
{
public:
    MongusFilter() : Filter()
    {
    }

    static void* create();
    static int32_t destroy(void*);
    std::string getName() const;

private:
    size_t m_rows;
    size_t m_cols;
    double m_cellSize;
    double m_k;
    int m_l;
    BOX2D m_bounds;
    DimRange m_ignored;
    bool m_lastOnly;
    std::string m_dir;

    virtual void addDimensions(PointLayoutPtr layout);
    virtual void addArgs(ProgramArgs& args);
    virtual void ready(PointTableRef table);
    size_t getColIndex(double x, double cell_size);
    size_t getRowIndex(double y, double cell_size);
    void writeControl(std::vector<double> cx, std::vector<double> cy,
                      std::vector<double> cz, std::string filename);
    void downsampleMin(std::vector<double>* cx, std::vector<double>* cy,
                       std::vector<double>* cz, std::vector<double>* dcx,
                       std::vector<double>* dcy, std::vector<double>* dcz,
                       double cell_size, size_t nr, size_t nc);
    std::vector<PointId> processGround(PointViewPtr view);
    virtual PointViewSet run(PointViewPtr view);

    MongusFilter& operator=(const MongusFilter&); // not implemented
    MongusFilter(const MongusFilter&);            // not implemented
};

} // namespace pdal
