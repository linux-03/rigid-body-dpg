#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include "rigid_body.h"
#include "rigid_body_system.h"
#include "Elements/Beam.h"
#include "Elements/Connector.h"
#include "Elements/Spring.h"
#include "equation_rbs.h"

namespace py = pybind11;

PYBIND11_MODULE(rigid_body_FEM, rbd) {

    // adds the BLA bindings as submodule, accessible as rigid_body.bla.Matrix etc.
    // https://github.com/pybind/pybind11/discussions/4027
    auto m = rbd.def_submodule("bla", "basic linear algebra");
    #include "bind_bla_obj.h"



    // the main bindings:
    rbd.doc() = "rigid body Shake and Rattle simulator";

    py::enum_<ConnectorType>(rbd, "ConnectorType", py::module_local())
      .value("FREE", ConnectorType::FREE)
      .value("FIX", ConnectorType::FIX)
      .value("SPHERICAL", ConnectorType::SPHERICAL)
      .export_values();

    py::class_<Connector>(rbd,"Connector", py::module_local())
      .def(py::init<>())
      .def(py::init<Vector<double>>())
      .def(py::init<Vector<double>, RigidBody&, ConnectorType>(),
           py::arg("pos"), py::arg("rb"), py::arg("type") = ConnectorType::FREE)
      .def_property("RefPos",[](Connector& c){return py::make_tuple(c.RefPosition(0),c.RefPosition(1),c.RefPosition(2));},
                             [](Connector&c, std::array<double,3> vals){c.RefPosition()(0) = vals[0]; c.RefPosition()(1) = vals[1]; c.RefPosition()(2) = vals[2];})
      .def_property("body_index", [](Connector& c){return c.BodyIndex();},
                                  [](Connector&c, size_t i){ c.BodyIndex() = i;})
      .def_property("type", [](Connector& c){return c.Type();},
                            [](Connector&c, ConnectorType t){ c.Type() = t;});

    py::class_<Beam>(rbd,"Beam", py::module_local())
      .def(py::init<>([](Connector a, Connector b){return Beam{a,b};}))
      .def_property_readonly("length", [](Beam& b){return b.Length();})
      .def_property_readonly("connectorA", [](Beam& b){return b.ConnectorA();})
      .def_property_readonly("connectorB",[](Beam& b){return b.ConnectorB();});

    py::class_<Spring>(rbd,"Spring", py::module_local())
      .def(py::init<>([](Connector a, Connector b){return Spring{a,b};})) // stiffness should be positive! (-k)
      .def_property_readonly("stiffness", [](Spring& b){return b.Stiffness();}) // stiffness should be positive! (-k)
      .def_property_readonly("connectorA", [](Spring& b){return b.ConnectorA();})
      .def_property_readonly("connectorB",[](Spring& b){return b.ConnectorB();});

    py::class_<RigidBody> (rbd, "RigidBody_FEM")
      .def(py::init<>())
      .def(py::init<Vector<double>, Vector<double>>())
      .def(py::init<double, Matrix<double>>())
      .def("asTuple",[](RigidBody& rb){
        // *column-major* transformation matrix as in https://threejs.org/docs/#api/en/math/Matrix4
        // converts Schöberl-style ordering of Q to column-major ordering of a three.js transformation matrix:
        Vector<double> q = rb.Vector_q();
        return py::make_tuple(q(3), q(6), q(9), 0, q(4), q(7), q(10), 0, q(5), q(8), q(11), 0, q(0), q(1), q(2), 1);
      })
      .def_property("q_trans",[](RigidBody& c){ return py::make_tuple(c.q_trans()(0), c.q_trans()(1), c.q_trans()(2));},
                              [](RigidBody&c, std::array<double,3> vals){c.q_trans() = {vals[0], vals[1], vals[2]};})
      .def_property("q", [](RigidBody& c){ 
                            Matrix<double> mat = c.q();
                            return mat;},
                         [](RigidBody&c, Matrix<double> mat){c.q() = mat;})
      .def_property("p_trans", [](RigidBody& c){ return py::make_tuple(c.p_trans()(0), c.p_trans()(1), c.p_trans()(2));},
                              [](RigidBody&c, std::array<double,3> vals){c.p_trans() = {vals[0], vals[1], vals[2]};})
      .def_property("p_skew", [](RigidBody& c){ return py::make_tuple(c.p_skew()(0), c.p_skew()(1), c.p_skew()(2));},
                              [](RigidBody&c, std::array<double,3> vals){c.p_skew() = {vals[0], vals[1], vals[2]};})
      .def_property_readonly("index", [](RigidBody& rb){ return rb.Index();})
      .def_property("vertices", [](RigidBody& rb){
        py::list py_list;
        for (const auto& v : rb.vertices()) {
          py_list.append(py::make_tuple(v[0], v[1], v[2]));
        }
        return py_list;
        }, 
        [](RigidBody& rb, py::list py_list){
        std::vector<std::array<double, 3>> vec;
        for (auto item : py_list) {
          py::tuple tup = py::cast<py::tuple>(item);
          vec.push_back({py::cast<double>(tup[0]), py::cast<double>(tup[1]), py::cast<double>(tup[2])});
        }
        rb.vertices_ = vec;
        })
      .def_property("normals", [](RigidBody& rb){
        py::list py_list;
        for (const auto& v : rb.normals()) {
          py_list.append(py::make_tuple(v[0], v[1], v[2]));
        }
        return py_list;
        }, 
        [](RigidBody& rb, py::list py_list){
        std::vector<std::array<double, 3>> vec;
        for (auto item : py_list) {
          py::tuple tup = py::cast<py::tuple>(item);
          vec.push_back({py::cast<double>(tup[0]), py::cast<double>(tup[1]), py::cast<double>(tup[2])});
        }
        rb.normals_ = vec;
        });
      //.def("set_normals_vertives")
          
    py::class_<RigidBodySystem> (rbd, "RBS_FEM")
      .def(py::init<>())
      .def("add", py::overload_cast<RigidBody&>(&RigidBodySystem::add))
      .def("add", py::overload_cast<Beam&>(&RigidBodySystem::add))
      .def("add", py::overload_cast<Spring&>(&RigidBodySystem::add))
      .def_property("gravity", [](RigidBodySystem& r){ auto& g = r.Gravity(); return py::make_tuple(g(0), g(1), g(2));},
                               [](RigidBodySystem& r, std::array<double,3> vals){ auto& g = r.Gravity(); g(0) = vals[0]; g(1) = vals[1]; g(2) = vals[2];})
      .def("saveState", &RigidBodySystem::SaveState)
      .def("expandState", &RigidBodySystem::ExpandState)
      .def("beams", py::overload_cast<>(&RigidBodySystem::Beams))
      .def("bodies", py::overload_cast<>(&RigidBodySystem::Bodies))
      .def("connectorPos", [](RigidBodySystem& r, Connector c){auto v = r.connectorPosition(c); return py::make_tuple(v(0),v(1),v(2));})
      .def("assemble", [](RigidBodySystem& r, double step_size){ return assemble(r, step_size);});

    py::class_<RigidBodySystemEquation> (rbd, "RigidBodyEquation")
      .def(py::init<RigidBodySystem&, double>())
      .def("step", &RigidBodySystemEquation::step);
    
    
    rbd.def("simulate",[](RigidBodySystem& rbs, double step_size, size_t steps) {simulate_rbs(rbs, step_size, steps);},
            py::arg("rbs"), py::arg("step_size"), py::arg("steps"));

    rbd.def("init_simulate", [](RigidBodySystem& rbs, double tend, double steps){});
}
