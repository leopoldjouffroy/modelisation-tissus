/*
**    TP CPE Lyon
**    Copyright (C) 2015 Damien Rohmer
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mesh_parametric_cloth.hpp"

#include "../lib/common/error_handling.hpp"
#include <cmath>

// debug
#include <iostream>

namespace cpe
{

void mesh_parametric_cloth::update_force()
{

    int const Nu = size_u();
    int const Nv = size_v();
    int const N_total = Nu*Nv;
    ASSERT_CPE(static_cast<int>(force_data.size()) == Nu*Nv , "Error of size");


    //Gravity
    static vec3 const g (0.0f,0.0f,-9.81f);
    vec3 const g_normalized = g/N_total;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            // Question 8
            //if(ku == 0 && kv == 0)
            //{
                //force(ku+1,kv) = vec3(0,0,0);
                //force(ku,kv+1) = vec3(0,0,0);!    
            //}
            //else
            //{
                //force(ku,kv) = g_normalized;
            //}

            force(ku,kv) = g_normalized;                
        }
    }

    //*************************************************************//
    // TO DO, Calculer les forces s'appliquant sur chaque sommet
    //*************************************************************//
    
    // Constantes de raideur
    float const K_structural = 0.1f;
    float const K_shear = 0.1f;
    float const K_bending = 0.1f;
    
    // Longueurs des ressorts au repos
    float const L01u = 1.0f/(Nu-1);
    float const L01v = 1.0f/(Nv-1);
    float const L02 = sqrt(pow(L01u,2) + pow(L01v,2));
    float const L03u = 2.0f*L01u;
    float const L03v = 2.0f*L01v;

    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {   
            /*** Structural ***/
            force(ku,kv) += force_voisin(K_structural, L01v, ku, kv, 0, 1);
            force(ku,kv) += force_voisin(K_structural, L01u, ku, kv, 1, 0);
            force(ku,kv) += force_voisin(K_structural, L01v, ku, kv, 0, -1);
            force(ku,kv) += force_voisin(K_structural, L01u, ku, kv, -1, 0);

            /*** Shearing ***/
            force(ku,kv) += force_voisin(K_shear, L02, ku, kv, 1, 1);
            force(ku,kv) += force_voisin(K_shear, L02, ku, kv, 1, -1);
            force(ku,kv) += force_voisin(K_shear, L02, ku, kv, -1, 1);
            force(ku,kv) += force_voisin(K_shear, L02, ku, kv, -1, -1);

            /*** Bending ***/
            force(ku,kv) += force_voisin(K_bending, L03v, ku, kv, 0, 2);
            force(ku,kv) += force_voisin(K_bending, L03u, ku, kv, 2, 0);
            force(ku,kv) += force_voisin(K_bending, L03v, ku, kv, 0, -2);
            force(ku,kv) += force_voisin(K_bending, L03u, ku, kv, -2, 0);
        }
    }

}

vec3 mesh_parametric_cloth::force_voisin(float K, float L0, int ku, int kv, int u_vois, int v_vois)
{
    int const Nu = size_u();
    int const Nv = size_v();

    int const ku_vois = ku+u_vois;
    int const kv_vois = kv+v_vois;

    vec3 f;

    if((ku_vois<0) || (ku_vois>Nu-1) || (kv_vois<0) || (kv_vois>Nv-1)){return f;}

    vec3 const u = vertex(ku,kv) - vertex(ku_vois,kv_vois);

    f = K * (L0-norm(u)) * u/norm(u);

    return f;
}

void mesh_parametric_cloth::integration_step(float const dt)
{
    ASSERT_CPE(speed_data.size() == force_data.size(),"Incorrect size");
    ASSERT_CPE(static_cast<int>(speed_data.size()) == size_vertex(),"Incorrect size");


    int const Nu = size_u();
    int const Nv = size_v();
    
    //*************************************************************//
    // TO DO: Calculer l'integration numerique des positions au cours de l'intervalle de temps dt.
    //*************************************************************//
    vec3 const centre = vec3(0.5f,0.05f,-1.1f);
    
    for(int ku=0; ku<Nu; ++ku)
    {
        for(int kv=0; kv<Nv; ++kv)
        {
            vec3& p = vertex(ku,kv);
            vec3& v = speed(ku,kv); 
            vec3& f = force(ku,kv);
            
            v = v + f*dt;
            p = p + v*dt;

            // interraction avec le plan
            if(p.z() < -1.101f && abs(p.x())<=1 && abs(p.y())<1)
            {
                p.z() = -1.100f;
                v.z() = 0;
                f.z() = 0;
            }

            // interraction avec la sphère
            if(norm(p-centre)<0.198f)
            {
                // vecteur directeur 
                vec3 vec_dir = (p-centre)/norm(p-centre);
                p = centre + 0.2f*vec_dir;
                v = vec3();
                f = vec3();
                //v = v + dot(v,vec_dir)*vec_dir; 
                //f = f + dot(f,vec_dir)*vec_dir;
            }
        }
    }

    //security check (throw exception if divergence is detected)
    static float const LIMIT=30.0f;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            vec3 const& p = vertex(ku,kv);

            if( norm(p) > LIMIT )
            {
                throw exception_divergence("Divergence of the system",EXCEPTION_PARAMETERS_CPE);
            }
        }
    }

}

void mesh_parametric_cloth::set_plane_xy_unit(int const size_u_param,int const size_v_param)
{
    mesh_parametric::set_plane_xy_unit(size_u_param,size_v_param);

    int const N = size_u()*size_v();
    speed_data.resize(N);
    force_data.resize(N);
}

vec3 const& mesh_parametric_cloth::speed(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3& mesh_parametric_cloth::speed(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3 const& mesh_parametric_cloth::force(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}

vec3& mesh_parametric_cloth::force(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}


}
