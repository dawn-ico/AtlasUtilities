//===--------------------------------------------------------------------------------*- C++ -*-===//
//                          _
//                         | |
//                       __| | __ ___      ___ ___
//                      / _` |/ _` \ \ /\ / / '_  |
//                     | (_| | (_| |\ V  V /| | | |
//                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include "toylibIO.h"

void debugDumpMesh(const toylib::Grid& mesh, const std::string prefix) {
  {
    char buf[256];
    sprintf(buf, "%sT.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(const auto& cellIt : mesh.faces()) {
      fprintf(fp, "%d %d %d\n", cellIt.vertex(0).id() + 1, cellIt.vertex(1).id() + 1,
              cellIt.vertex(2).id() + 1);
    }
    fclose(fp);
  }
  {
    char buf[256];
    sprintf(buf, "%sP.txt", prefix.c_str());
    FILE* fp = fopen(buf, "w+");
    for(const auto& nodeIt : mesh.vertices()) {
      fprintf(fp, "%f %f \n", nodeIt.x(), nodeIt.y());
    }
    fclose(fp);
  }
}

void dumpMesh(const toylib::Grid& m, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& e : m.edges()) {
    fprintf(fp, "%f %f %f %f\n", e.get().vertex(0).x(), e.get().vertex(0).y(),
            e.get().vertex(1).x(), e.get().vertex(1).y());
  }
  fclose(fp);
}

void dumpDualMesh(const toylib::Grid& m, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& e : m.edges()) {
    if(e.get().faces().size() != 2) {
      continue;
    }

    // This is WRONG!, leads to a dual mesh which is not orthogonal to
    // primal mesh
    // auto [xm1, ym1] = CellMidPoint(e.get().face(0));
    // auto [xm2, ym2] = CellMidPoint(e.get().face(1));

    auto [xm1, ym1] = CellCircumcenter(e.get().face(0));
    auto [xm2, ym2] = CellCircumcenter(e.get().face(1));
    fprintf(fp, "%f %f %f %f\n", xm1, ym1, xm2, ym2);
  }
  fclose(fp);
}

void dumpSparseData(const toylib::Grid& mesh, const toylib::SparseVertexData<double>& sparseData,
                    int level, int edgesPerVertex, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& v : mesh.vertices()) {
    int sparse_idx = 0;
    for(const auto& e : v.edges()) {
      double val = sparseData(v, sparse_idx, level);
      auto [emx, emy] = EdgeMidpoint(*e);
      double dx = emx - v.x();
      double dy = emy - v.y();
      fprintf(fp, "%f %f %f\n", v.x() + 0.5 * dx, v.y() + 0.5 * dy, val);
      sparse_idx++;
    }
  }
}

void dumpSparseData(const toylib::Grid& mesh, const toylib::SparseFaceData<double>& sparseData,
                    int level, int edgesPerCell, const std::string& fname) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(const auto& c : mesh.faces()) {
    int sparse_idx = 0;
    auto [cmx, cmy] = CellCircumcenter(c);
    for(const auto& e : c.edges()) {
      double val = sparseData(c, sparse_idx, level);
      auto [emx, emy] = EdgeMidpoint(*e);
      double dx = emx - cmx;
      double dy = emy - cmy;
      fprintf(fp, "%f %f %f\n", cmx + 0.5 * dx, cmy + 0.5 * dy, val);
      sparse_idx++;
    }
  }
}

void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::EdgeData<double>& field, int level,
               std::optional<toylib::edge_color> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& e : mesh.edges()) {
    if(color.has_value() && e.get().color() != color.value()) {
      continue;
    }
    auto [x, y] = EdgeMidpoint(e);
    fprintf(fp, "%f %f %f\n", x, y, std::isfinite(field(e, level)) ? field(e, level) : 0.);
  }
  fclose(fp);
}

void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::EdgeData<double>& field_x, const toylib::EdgeData<double>& field_y,
               int level, std::optional<toylib::edge_color> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& e : mesh.edges()) {
    if(color.has_value() && e.get().color() != color.value()) {
      continue;
    }
    auto [x, y] = EdgeMidpoint(e);
    fprintf(fp, "%f %f %f %f\n", x, y, field_x(e, level), field_y(e, level));
  }
  fclose(fp);
}

void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::FaceData<double>& field, int level,
               std::optional<toylib::face_color> color) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& c : mesh.faces()) {
    if(color.has_value() && c.color() != color.value()) {
      continue;
    }
    auto [x, y] = CellCircumcenter(c);
    fprintf(fp, "%f %f %f\n", x, y, field(c, level));
  }
  fclose(fp);
}

void dumpField(const std::string& fname, const toylib::Grid& mesh,
               const toylib::VertexData<double>& field, int level) {
  FILE* fp = fopen(fname.c_str(), "w+");
  for(auto& v : mesh.vertices()) {
    fprintf(fp, "%f %f %f\n", v.x(), v.y(), field(v, level));
  }
  fclose(fp);
}