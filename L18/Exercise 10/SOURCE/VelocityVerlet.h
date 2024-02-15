void PBC_no_move(int idp){
    parts[idp].x = P_Img(parts[idp].ox, mySys.box_x);
    parts[idp].y = P_Img(parts[idp].oy, mySys.box_y);
    parts[idp].z = P_Img(parts[idp].oz, mySys.box_z);
}

void set_forces(int idp)
{
    parts[idp].fx = 0.;
    parts[idp].fy = 0.;
    parts[idp].fz = 0.;
}

double compute_force(double r){
    return (r < mySys.sigma_cut) ? 24.*mySys.eps*(2.*pow(mySys.sigma/r, 12.) - pow(mySys.sigma/r, 6.))/pow(r, 2.) : 0.;
}

void Velocity_Verlet(double dt){
    
    double r_ij = 0;
    for(int i = 0; i < mySys.NPart; i++){

        // momentum translation by an amount (half a kick)
        parts[i].ovx += dt* parts[i].fx/2.;
        parts[i].ovy += dt* parts[i].fy/2.;
        parts[i].ovz += dt* parts[i].fz/2.;

        // position translation by an amount (drift)
        parts[i].ox += parts[i].ovx * dt;
        parts[i].oy += parts[i].ovy * dt;
        parts[i].oz += parts[i].ovz * dt;
        
    }
    
    // periodic boundary conditions: confine particles inside the box
    // update positions
    for(int l = 0; l < mySys.NPart; l++){
        PBC_no_move(l);
        set_forces(l); // set forces to 0
    }  
    
    // forces update:
    double dxu = 0, dyu = 0, dzu = 0;
    double force = 0.;
    for(int i = 0; i < mySys.NPart; i++){
        for(int j = 0; j < mySys.NPart; j++){
            if (j == i) continue;    
            // compute new forces for the updated positions:
            r_ij = distance_updated(parts[i], parts[j]); 
            dxu = parts[i].x - parts[j].x;
            dyu = parts[i].y - parts[j].y;
            dzu = parts[i].z - parts[j].z;
            
            // Periodic Boundary Conditions
            dxu = dxu - round(dxu/mySys.box_x)*mySys.box_x;
            dyu = dyu - round(dyu/mySys.box_y)*mySys.box_y;
            dzu = dzu - round(dzu/mySys.box_z)*mySys.box_z;

            force = compute_force(r_ij);
            parts[i].fx += force*dxu; 
            parts[i].fy += force*dyu; 
            parts[i].fz += force*dzu; 
        }

        // momentum translation by an amount (half a kick)
        parts[i].ovx += dt* parts[i].fx/2.;
        parts[i].ovy += dt* parts[i].fy/2.;
        parts[i].ovz += dt* parts[i].fz/2.;
     
        // update velocities (momenta)
        parts[i].vx = parts[i].ovx; parts[i].vy = parts[i].ovy; parts[i].vz = parts[i].ovz;
    }
}
