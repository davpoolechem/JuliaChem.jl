#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>

#include <unistd.h>
#include <iostream>
#include <vector>

typedef int64_t julia_int;

//--------------------------------//
//-- Map libint2::Atom to Julia --// 
//--------------------------------//
template<> struct jlcxx::IsMirroredType<libint2::Atom> : std::false_type { };

julia_int get_atomic_number(const libint2::Atom& atom) { 
  return atom.atomic_number; 
}
julia_int set_atomic_number(libint2::Atom& atom, julia_int new_an) { 
  atom.atomic_number = new_an;
  return get_atomic_number(atom); 
}

double get_x(const libint2::Atom& atom) { return atom.x; }
double set_x(libint2::Atom& atom, double new_x) {
  atom.x = new_x; 
  return get_x(atom); 
}

double get_y(const libint2::Atom& atom) { return atom.y; }
double set_y(libint2::Atom& atom, double new_y) {
  atom.y = new_y; 
  return get_y(atom); 
}

double get_z(const libint2::Atom& atom) { return atom.z; }
double set_z(libint2::Atom& atom, double new_z) {
  atom.z = new_z; 
  return get_z(atom); 
}

//------------------------------------//
//-- Map libint2::BasisSet to Julia --// 
//------------------------------------//
template<> struct jlcxx::IsMirroredType<std::vector<libint2::Atom> > : std::false_type { };
template<> struct jlcxx::IsMirroredType<libint2::BasisSet> : std::false_type { };

//------------------------------------------------------------------------//
//-- C++ JERI engine: small wrapper allowing LibInt to be used in Julia --//
//------------------------------------------------------------------------//
class Engine {
  libint2::BasisSet m_basis_set;
  
  libint2::Engine m_kinetic_eng;
  libint2::Engine m_overlap_eng;
  libint2::Engine m_nuc_attr_eng;

public:
  //-- ctors and dtors --//
  Engine() { initialize(); };
  Engine(const std::vector<libint2::Atom>& t_atoms, const std::string& basis) 
    : m_basis_set(basis, t_atoms) 
  { 
    initialize();
 
    m_overlap_eng = libint2::Engine(libint2::Operator::overlap, 
      m_basis_set.max_nprim(), m_basis_set.max_l(), 0);
    //m_overlap_eng.set_params();

    m_kinetic_eng = libint2::Engine(libint2::Operator::kinetic, 
      m_basis_set.max_nprim(), m_basis_set.max_l(), 0);
    //m_kinetic_eng.set_params();

    m_nuc_attr_eng = libint2::Engine(libint2::Operator::nuclear, 
      m_basis_set.max_nprim(), m_basis_set.max_l(), 0);
    m_nuc_attr_eng.set_params(libint2::make_point_charges(t_atoms));
  }

  ~Engine() { finalize(); };

  //-- setters and getters --//
  libint2::BasisSet basis() { return m_basis_set; }
  //libint2::Engine engine() { return m_engine; }

  //-- member functions --//
  void compute_overlap_block(jlcxx::ArrayRef<double> S_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_overlap_eng.compute(m_basis_set[ash-1], m_basis_set[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      S_block[i] = m_overlap_eng.results()[0][i];
    }
  }
 
  void compute_kinetic_block(jlcxx::ArrayRef<double> T_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_kinetic_eng.compute(m_basis_set[ash-1], m_basis_set[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      T_block[i] = m_kinetic_eng.results()[0][i];
    }
  }

  void compute_nuc_attr_block(jlcxx::ArrayRef<double> V_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_nuc_attr_eng.compute(m_basis_set[ash-1], m_basis_set[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      V_block[i] = m_nuc_attr_eng.results()[0][i];
    }
  }

private:
  //-- private member functions --//
  void initialize() {
    libint2::initialize(); 
  }

  void finalize() {
    libint2::finalize();
  }
};

JLCXX_MODULE define_jeri(jlcxx::Module& mod) {
  //-- atom information --//
  mod.add_type<libint2::Atom>("Atom");
  jlcxx::stl::apply_stl<libint2::Atom>(mod);

  mod.method("atomic_number", &get_atomic_number);
  mod.method("atomic_number", &set_atomic_number);

  mod.method("x",&get_x);
  mod.method("x",&set_x);

  mod.method("y",&get_y);
  mod.method("y",&set_y);

  mod.method("z",&get_z);
  mod.method("z",&set_z);

  //-- basis set information --//
  //mod.add_type<std::vector<libint2::Atom> >("Atoms");
  mod.add_type<libint2::BasisSet>("BasisSet");

  //-- engine information --//
  mod.add_type<libint2::Engine>("LibIntEngine");
  mod.add_type<Engine>("Engine")
    .constructor<const std::vector<libint2::Atom>&, const std::string& >()
    .method("basis", &Engine::basis)
    .method("compute_overlap_block", &Engine::compute_overlap_block)
    .method("compute_kinetic_block", &Engine::compute_kinetic_block)
    .method("compute_nuc_attr_block", &Engine::compute_nuc_attr_block);
} 
