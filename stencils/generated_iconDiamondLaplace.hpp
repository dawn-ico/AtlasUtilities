//---- Preprocessor defines ----
#define DAWN_GENERATED 1
#undef DAWN_BACKEND_T
#define DAWN_BACKEND_T CXXNAIVEICO
#include "interfaces/unstructured_interface.hpp"

//---- Globals ----

//---- Stencils ----
namespace dawn_generated {
namespace cxxnaiveico {
template <typename LibTag>
class ICON_laplacian_diamond_stencil {
private:
  struct stencil_175 {
    dawn::mesh_t<LibTag> const& m_mesh;
    int m_k_size;
    dawn::edge_field_t<LibTag, double>& m_diff_multfac_smag;
    dawn::edge_field_t<LibTag, double>& m_tangent_orientation;
    dawn::edge_field_t<LibTag, double>& m_inv_primal_edge_length;
    dawn::edge_field_t<LibTag, double>& m_inv_vert_vert_length;
    dawn::vertex_field_t<LibTag, double>& m_u_vert;
    dawn::vertex_field_t<LibTag, double>& m_v_vert;
    dawn::sparse_edge_field_t<LibTag, double>& m_primal_normal_x;
    dawn::sparse_edge_field_t<LibTag, double>& m_primal_normal_y;
    dawn::sparse_edge_field_t<LibTag, double>& m_dual_normal_x;
    dawn::sparse_edge_field_t<LibTag, double>& m_dual_normal_y;
    dawn::sparse_edge_field_t<LibTag, double>& m_vn_vert;
    dawn::edge_field_t<LibTag, double>& m_vn;
    dawn::edge_field_t<LibTag, double>& m_dvt_tang;
    dawn::edge_field_t<LibTag, double>& m_dvt_norm;
    dawn::edge_field_t<LibTag, double>& m_kh_smag_1;
    dawn::edge_field_t<LibTag, double>& m_kh_smag_2;
    dawn::edge_field_t<LibTag, double>& m_kh_smag;
    dawn::edge_field_t<LibTag, double>& m_nabla2;

  public:
    stencil_175(
        dawn::mesh_t<LibTag> const& mesh, int k_size,
        dawn::edge_field_t<LibTag, double>& diff_multfac_smag,
        dawn::edge_field_t<LibTag, double>& tangent_orientation,
        dawn::edge_field_t<LibTag, double>& inv_primal_edge_length,
        dawn::edge_field_t<LibTag, double>& inv_vert_vert_length,
        dawn::vertex_field_t<LibTag, double>& u_vert, dawn::vertex_field_t<LibTag, double>& v_vert,
        dawn::sparse_edge_field_t<LibTag, double>& primal_normal_x,
        dawn::sparse_edge_field_t<LibTag, double>& primal_normal_y,
        dawn::sparse_edge_field_t<LibTag, double>& dual_normal_x,
        dawn::sparse_edge_field_t<LibTag, double>& dual_normal_y,
        dawn::sparse_edge_field_t<LibTag, double>& vn_vert, dawn::edge_field_t<LibTag, double>& vn,
        dawn::edge_field_t<LibTag, double>& dvt_tang, dawn::edge_field_t<LibTag, double>& dvt_norm,
        dawn::edge_field_t<LibTag, double>& kh_smag_1,
        dawn::edge_field_t<LibTag, double>& kh_smag_2, dawn::edge_field_t<LibTag, double>& kh_smag,
        dawn::edge_field_t<LibTag, double>& nabla2)
        : m_mesh(mesh), m_k_size(k_size), m_diff_multfac_smag(diff_multfac_smag),
          m_tangent_orientation(tangent_orientation),
          m_inv_primal_edge_length(inv_primal_edge_length),
          m_inv_vert_vert_length(inv_vert_vert_length), m_u_vert(u_vert), m_v_vert(v_vert),
          m_primal_normal_x(primal_normal_x), m_primal_normal_y(primal_normal_y),
          m_dual_normal_x(dual_normal_x), m_dual_normal_y(dual_normal_y), m_vn_vert(vn_vert),
          m_vn(vn), m_dvt_tang(dvt_tang), m_dvt_norm(dvt_norm), m_kh_smag_1(kh_smag_1),
          m_kh_smag_2(kh_smag_2), m_kh_smag(kh_smag), m_nabla2(nabla2) {}

    ~stencil_175() {}

    void sync_storages() {}

    void run() {
      using dawn::deref;
      {
        for(int k = 0 + 0; k <= (m_k_size == 0 ? 0 : (m_k_size - 1)) + 0 + 0; ++k) {
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            {
              int for_loop_idx = 0;
              for(auto inner_loc :
                  getNeighbors(LibTag{}, m_mesh,
                               std::vector<dawn::LocationType>{dawn::LocationType::Edges,
                                                               dawn::LocationType::Cells,
                                                               dawn::LocationType::Vertices},
                               loc)) {
                m_vn_vert(deref(LibTag{}, loc), for_loop_idx, k + 0) =
                    ((m_u_vert(deref(LibTag{}, inner_loc), k + 0) *
                      m_primal_normal_x(deref(LibTag{}, loc), for_loop_idx, k + 0)) +
                     (m_v_vert(deref(LibTag{}, inner_loc), k + 0) *
                      m_primal_normal_y(deref(LibTag{}, loc), for_loop_idx, k + 0)));
                for_loop_idx++;
              }
            }
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            {
              int sparse_dimension_idx0 = 0;
              m_dvt_tang(deref(LibTag{}, loc), k + 0) = reduce(
                  LibTag{}, m_mesh, loc, (::dawn::float_type)0.0,
                  std::vector<dawn::LocationType>{dawn::LocationType::Edges,
                                                  dawn::LocationType::Cells,
                                                  dawn::LocationType::Vertices},
                  [&](auto& lhs, auto red_loc1, auto const& weight) {
                    lhs += weight *
                           ((m_u_vert(deref(LibTag{}, red_loc1), k + 0) *
                             m_dual_normal_x(deref(LibTag{}, loc), sparse_dimension_idx0, k + 0)) +
                            (m_v_vert(deref(LibTag{}, red_loc1), k + 0) *
                             m_dual_normal_y(deref(LibTag{}, loc), sparse_dimension_idx0, k + 0)));
                    sparse_dimension_idx0++;
                    return lhs;
                  },
                  std::vector<::dawn::float_type>({(::dawn::float_type)-1.0,
                                                   (::dawn::float_type)1.0, (::dawn::float_type)0.0,
                                                   (::dawn::float_type)0.0}));
            }
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            m_dvt_tang(deref(LibTag{}, loc), k + 0) =
                (m_dvt_tang(deref(LibTag{}, loc), k + 0) *
                 m_tangent_orientation(deref(LibTag{}, loc), k + 0));
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            {
              int sparse_dimension_idx0 = 0;
              m_dvt_norm(deref(LibTag{}, loc), k + 0) = reduce(
                  LibTag{}, m_mesh, loc, (::dawn::float_type)0.0,
                  std::vector<dawn::LocationType>{dawn::LocationType::Edges,
                                                  dawn::LocationType::Cells,
                                                  dawn::LocationType::Vertices},
                  [&](auto& lhs, auto red_loc1, auto const& weight) {
                    lhs += weight *
                           ((m_u_vert(deref(LibTag{}, red_loc1), k + 0) *
                             m_dual_normal_x(deref(LibTag{}, loc), sparse_dimension_idx0, k + 0)) +
                            (m_v_vert(deref(LibTag{}, red_loc1), k + 0) *
                             m_dual_normal_y(deref(LibTag{}, loc), sparse_dimension_idx0, k + 0)));
                    sparse_dimension_idx0++;
                    return lhs;
                  },
                  std::vector<::dawn::float_type>({(::dawn::float_type)0.0, (::dawn::float_type)0.0,
                                                   (::dawn::float_type)-1.0,
                                                   (::dawn::float_type)1.0}));
            }
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            {
              int sparse_dimension_idx0 = 0;
              m_kh_smag_1(deref(LibTag{}, loc), k + 0) = reduce(
                  LibTag{}, m_mesh, loc, (::dawn::float_type)0.0,
                  std::vector<dawn::LocationType>{dawn::LocationType::Edges,
                                                  dawn::LocationType::Cells,
                                                  dawn::LocationType::Vertices},
                  [&](auto& lhs, auto red_loc1, auto const& weight) {
                    lhs += weight * m_vn_vert(deref(LibTag{}, loc), sparse_dimension_idx0, k + 0);
                    sparse_dimension_idx0++;
                    return lhs;
                  },
                  std::vector<::dawn::float_type>({(::dawn::float_type)-1.0,
                                                   (::dawn::float_type)1.0, (::dawn::float_type)0.0,
                                                   (::dawn::float_type)0.0}));
            }
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            m_kh_smag_1(deref(LibTag{}, loc), k + 0) =
                (((m_kh_smag_1(deref(LibTag{}, loc), k + 0) *
                   m_tangent_orientation(deref(LibTag{}, loc), k + 0)) *
                  m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0)) +
                 (m_dvt_norm(deref(LibTag{}, loc), k + 0) *
                  m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0)));
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            m_kh_smag_1(deref(LibTag{}, loc), k + 0) = (m_kh_smag_1(deref(LibTag{}, loc), k + 0) *
                                                        m_kh_smag_1(deref(LibTag{}, loc), k + 0));
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            {
              int sparse_dimension_idx0 = 0;
              m_kh_smag_2(deref(LibTag{}, loc), k + 0) = reduce(
                  LibTag{}, m_mesh, loc, (::dawn::float_type)0.0,
                  std::vector<dawn::LocationType>{dawn::LocationType::Edges,
                                                  dawn::LocationType::Cells,
                                                  dawn::LocationType::Vertices},
                  [&](auto& lhs, auto red_loc1, auto const& weight) {
                    lhs += weight * m_vn_vert(deref(LibTag{}, loc), sparse_dimension_idx0, k + 0);
                    sparse_dimension_idx0++;
                    return lhs;
                  },
                  std::vector<::dawn::float_type>({(::dawn::float_type)0.0, (::dawn::float_type)0.0,
                                                   (::dawn::float_type)-1.0,
                                                   (::dawn::float_type)1.0}));
            }
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            m_kh_smag_2(deref(LibTag{}, loc), k + 0) =
                ((m_kh_smag_2(deref(LibTag{}, loc), k + 0) *
                  m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0)) -
                 (m_dvt_tang(deref(LibTag{}, loc), k + 0) *
                  m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0)));
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            m_kh_smag_2(deref(LibTag{}, loc), k + 0) = (m_kh_smag_2(deref(LibTag{}, loc), k + 0) *
                                                        m_kh_smag_2(deref(LibTag{}, loc), k + 0));
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            m_kh_smag(deref(LibTag{}, loc), k + 0) =
                (m_diff_multfac_smag(deref(LibTag{}, loc), k + 0) *
                 (m_kh_smag_1(deref(LibTag{}, loc), k + 0) +
                  m_kh_smag_2(deref(LibTag{}, loc), k + 0)));
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            {
              int sparse_dimension_idx0 = 0;
              m_nabla2(deref(LibTag{}, loc), k + 0) = reduce(
                  LibTag{}, m_mesh, loc, (::dawn::float_type)0.0,
                  std::vector<dawn::LocationType>{dawn::LocationType::Edges,
                                                  dawn::LocationType::Cells,
                                                  dawn::LocationType::Vertices},
                  [&](auto& lhs, auto red_loc1, auto const& weight) {
                    lhs += weight * ((::dawn::float_type)4.0 *
                                     m_vn_vert(deref(LibTag{}, loc), sparse_dimension_idx0, k + 0));
                    sparse_dimension_idx0++;
                    return lhs;
                  },
                  std::vector<::dawn::float_type>(
                      {(m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0) *
                        m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0)),
                       (m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0) *
                        m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0)),
                       (m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0) *
                        m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0)),
                       (m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0) *
                        m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0))}));
            }
          }
          for(auto const& loc : getEdges(LibTag{}, m_mesh)) {
            m_nabla2(deref(LibTag{}, loc), k + 0) =
                (m_nabla2(deref(LibTag{}, loc), k + 0) -
                 ((((::dawn::float_type)8.0 * m_vn(deref(LibTag{}, loc), k + 0)) *
                   (m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0) *
                    m_inv_primal_edge_length(deref(LibTag{}, loc), k + 0))) +
                  (((::dawn::float_type)8.0 * m_vn(deref(LibTag{}, loc), k + 0)) *
                   (m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0) *
                    m_inv_vert_vert_length(deref(LibTag{}, loc), k + 0)))));
          }
        }
      }
      sync_storages();
    }
  };
  static constexpr const char* s_name = "ICON_laplacian_diamond_stencil";
  stencil_175 m_stencil_175;

public:
  ICON_laplacian_diamond_stencil(const ICON_laplacian_diamond_stencil&) = delete;

  // Members

  ICON_laplacian_diamond_stencil(
      const dawn::mesh_t<LibTag>& mesh, int k_size,
      dawn::edge_field_t<LibTag, double>& diff_multfac_smag,
      dawn::edge_field_t<LibTag, double>& tangent_orientation,
      dawn::edge_field_t<LibTag, double>& inv_primal_edge_length,
      dawn::edge_field_t<LibTag, double>& inv_vert_vert_length,
      dawn::vertex_field_t<LibTag, double>& u_vert, dawn::vertex_field_t<LibTag, double>& v_vert,
      dawn::sparse_edge_field_t<LibTag, double>& primal_normal_x,
      dawn::sparse_edge_field_t<LibTag, double>& primal_normal_y,
      dawn::sparse_edge_field_t<LibTag, double>& dual_normal_x,
      dawn::sparse_edge_field_t<LibTag, double>& dual_normal_y,
      dawn::sparse_edge_field_t<LibTag, double>& vn_vert, dawn::edge_field_t<LibTag, double>& vn,
      dawn::edge_field_t<LibTag, double>& dvt_tang, dawn::edge_field_t<LibTag, double>& dvt_norm,
      dawn::edge_field_t<LibTag, double>& kh_smag_1, dawn::edge_field_t<LibTag, double>& kh_smag_2,
      dawn::edge_field_t<LibTag, double>& kh_smag, dawn::edge_field_t<LibTag, double>& nabla2)
      : m_stencil_175(mesh, k_size, diff_multfac_smag, tangent_orientation, inv_primal_edge_length,
                      inv_vert_vert_length, u_vert, v_vert, primal_normal_x, primal_normal_y,
                      dual_normal_x, dual_normal_y, vn_vert, vn, dvt_tang, dvt_norm, kh_smag_1,
                      kh_smag_2, kh_smag, nabla2) {}

  void run() {
    m_stencil_175.run();
    ;
  }
};
} // namespace cxxnaiveico
} // namespace dawn_generated
