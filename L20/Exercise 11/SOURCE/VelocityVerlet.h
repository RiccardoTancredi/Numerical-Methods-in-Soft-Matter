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

void set_old_forces(int idp)
{
    parts[idp].ofx = parts[idp].fx;
    parts[idp].ofy = parts[idp].fy;
    parts[idp].ofz = parts[idp].fz;
}

// double compute_force(double r){
//     return (r < mySys.sigma_cut) ? 24.*mySys.eps*(2.*pow(mySys.sigma/r, 12.) - pow(mySys.sigma/r, 6.))/pow(r, 2.) : 0.;
// }

double compute_force(){
    // Work in 2D here
    return -mySys.K;
}



void overdamped_Langevin_equation(double dt){
    
    double x1 = 0, x2 = 0, x3 = 0;  // Box-Muller variables
    for(int i = 0; i < mySys.NPart; i++){

        // position translation by an amount (drift)
        BoxMuller(&x1, &x2); BoxMuller(&x3, &x2);
        parts[i].ox += parts[i].fx * dt/(mySys.m*mySys.gamma) + sqrt((2.*mySys.T*dt)/(mySys.m*mySys.gamma))*x1;
        parts[i].oy += parts[i].fy * dt/(mySys.m*mySys.gamma) + sqrt((2.*mySys.T*dt)/(mySys.m*mySys.gamma))*x2;
        if (mySys.exercise != 2) parts[i].oz += parts[i].fz * dt/(mySys.m*mySys.gamma) + sqrt((2.*mySys.T*dt)/(mySys.m*mySys.gamma))*x3;
    }
    
    // periodic boundary conditions: confine particles inside the box
    // update positions
    for(int l = 0; l < mySys.NPart; l++){
        PBC_no_move(l);
        set_forces(l); // set forces to 0
    }  
}

void Velocity_Verlet_underdamped(double dt){
    
    // double r_ij = 0;
    double x1 = 0, x2 = 0, x3 = 0;  // Box-Muller variables
    double theta1 = 0, theta2 = 0, theta3 = 0;  // Box-Muller variables

    // update C(t)
    for(int i = 0; i < mySys.NPart; i++){
        BoxMuller(&x1, &x2); BoxMuller(&x3, &x2);
        BoxMuller(&theta1, &theta2); BoxMuller(&theta3, &theta2);
        parts[i].Cx = dt*dt/2. * (parts[i].fx - mySys.gamma*parts[i].ovx) + sqrt((2.*mySys.T*mySys.gamma)/(mySys.m))*pow(dt, 3./2.)*(0.5*x1+theta1/(2.*sqrt(3)));
        parts[i].Cy = dt*dt/2. * (parts[i].fy - mySys.gamma*parts[i].ovy) + sqrt((2.*mySys.T*mySys.gamma)/(mySys.m))*pow(dt, 3./2.)*(0.5*x2+theta2/(2.*sqrt(3)));
        if (mySys.exercise != 2) parts[i].Cz = dt*dt/2. * (parts[i].fz - mySys.gamma*parts[i].ovz) + sqrt((2.*mySys.T*mySys.gamma)/(mySys.m))*pow(dt, 3./2.)*(0.5*x3+theta3/(2.*sqrt(3)));
    }

    // update postions
    for(int i = 0; i < mySys.NPart; i++){

        // position translation
        parts[i].ox += parts[i].ovx * dt + parts[i].Cx;
        parts[i].oy += parts[i].ovy * dt + parts[i].Cy;
        if (mySys.exercise != 2) parts[i].oz += parts[i].ovz * dt + parts[i].Cz;
    }
    
    // periodic boundary conditions: confine particles inside the box
    // update positions
    for(int l = 0; l < mySys.NPart; l++){
        PBC_no_move(l);
        set_old_forces(l); // set forces to old forces
        set_forces(l); // set forces to 0
    }  
    
    // forces update:
    double dxu = 0, dyu = 0, dzu = 0;
    double force = 0.;

    for(int i = 0; i < mySys.NPart; i++){
        if (mySys.exercise == 2){
            // compute new forces for the updated positions:
            dxu = parts[i].x;
            dyu = parts[i].y;
            
            // Periodic Boundary Conditions
            dxu = dxu - round(dxu/mySys.box_x)*mySys.box_x;
            dyu = dyu - round(dyu/mySys.box_y)*mySys.box_y;

            force = compute_force();
            parts[i].fx += force*dxu; 
            parts[i].fy += force*dyu; 

            if (mySys.spring_2_x != 0 || mySys.spring_2_y != 0){
                // There is another spring
                dxu = parts[i].ox - mySys.spring_2_x;
                dyu = parts[i].oy - mySys.spring_2_y;
                // dzu = parts[i].oz; 

                // Periodic Boundary Conditions
                dxu = dxu - round(dxu/mySys.box_x)*mySys.box_x;
                dyu = dyu - round(dyu/mySys.box_y)*mySys.box_y;

                parts[i].fx += force*dxu; 
                parts[i].fy += force*dyu;
            }

        }

        BoxMuller(&x1, &x2); BoxMuller(&x3, &x2);

        // momentum translation by an amount (half a kick)
        parts[i].ovx += dt* (parts[i].fx + parts[i].ofx)/2. - dt*mySys.gamma*parts[i].ovx + sqrt((2.*mySys.T*mySys.gamma*dt)/(mySys.m))*x1 - mySys.gamma*parts[i].Cx;
        parts[i].ovy += dt* (parts[i].fy + parts[i].ofy)/2. - dt*mySys.gamma*parts[i].ovy + sqrt((2.*mySys.T*mySys.gamma*dt)/(mySys.m))*x2 - mySys.gamma*parts[i].Cy;
        if (mySys.exercise != 2) parts[i].ovz += dt* (parts[i].fz + parts[i].ofz)/2. - dt*mySys.gamma*parts[i].ovz + sqrt((2.*mySys.T*mySys.gamma*dt)/(mySys.m))*x3 - mySys.gamma*parts[i].Cz;
     
        // update velocities (momenta)
        parts[i].vx = parts[i].ovx; parts[i].vy = parts[i].ovy; 
        if (mySys.exercise != 2) parts[i].vz = parts[i].ovz;
    }
}
