package comp559.fluid;

import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2f;
import javax.vecmath.Vector2f;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;

/**
 * Eulerian fluid simulation class. 
 * 
 * This follows closely the solver documented in J. Stam GDC 2003, specifically
 * this uses a non staggered grid, with an extra set of grid cells to help enforce
 * boundary conditions
 * 
 * @author kry
 */
public class Fluid {
    
    final static int DIM = 2;
    
    /** velocity, non-staggered, packed (first index is the dimension), public for achievement checking */ 
    public float[][] U0;
    
    /** temporary velocity variable */
    private float[][] U1;          
       
    /** temperature (packed)*/
    public float[] temperature0;
    
    /** temperature (packed) temporary variable*/
    private float[] temperature1;
       
    private float[] tempTemp;
    private float[][] tempU;
    private IntParameter Nval = new IntParameter( "grid size (requires reset)", 16, 1, 256 );
    
    /** Number of grid cells (not counting extra boundary cells */
    public int N = 16;
    
    /** Dimension of each grid cell */
    public float dx = 1;
    
    /** time elapsed in the fluid simulation */
    public double elapsed;
    
    private float[] p;
    private float[] div;
    public float[] diagonalA;
    public float[] plusiA;
    public float[] plusjA;
    public float[] preconditioner;
    public float[] m;
    public float[] y;
    public float[] x;
    public float[] r;
    public float[] z;
    public float[] s;
    public float dt;
    /**
     * Sources of heat and cold
     */
    
    public List<Source> sources = new LinkedList<Source>();
    
    /**
     * initialize memory
     */
    public void setup() {
        elapsed = 0;
        N = Nval.getValue();        
        dx = 1.0f / N; // we choose the domain size here to be 1 unit square!
        
        int np2s = (N+2)*(N+2);
        U0 = new float[2][np2s];
        U1 = new float[2][np2s];
//        for (int i=0;i<np2s;i++){
//        	U0[0][i]=0;
//        }
        temperature0 = new float[np2s];
        temperature1 = new float[np2s];
        
        //  you'll also want to implement work arrays for the divergence and the result of the Poisson solve

        p = new float[np2s];
        div = new float[np2s];
        diagonalA = new float[np2s];
        plusiA = new float[np2s];
        plusjA = new float[np2s];
        preconditioner = new float[np2s];
        m = new float[np2s];
        y = new float[np2s];
        x = new float[np2s];
        r = new float[np2s];
        z = new float[np2s];
        s = new float[np2s];
    }

    /**
     * Compute the index 
     * @param i 
     * @param j 
     * @return the index in the flattened array
     */
    public int IX( int i, int j ) {
        return i*(N+2) + j;
    }
    
    /** 
     * Gets the velocity at the given point using interpolation 
     * @param x
     * @param vel
     */
    public void getVelocity( Tuple2f x, Tuple2f vel ) {
        getVelocity( x, U0, vel );
    }
    
    /** 
     * Gets the velocity in the provided velocity field at the given point using interpolation
     * @param x
     * @param U
     * @param vel
     */
    private void getVelocity( Tuple2f x, float[][] U, Tuple2f vel ) {
        vel.x = interpolateU( x, U[0] );
        vel.y = interpolateV( x, U[1] );
    }
    
    /**
     * Interpolates the given scalar field
     * @param x
     * @param s
     * @return interpolated value
     */
    public float interpolate( Tuple2f x, float[] s ) {
        int i,j;
        float wx, wy;
        i = (int) Math.floor( x.x / dx - 0.5 );
        j = (int) Math.floor( x.y / dx - 0.5 );
        wx = x.x / dx - 0.5f - i;
        wy = x.y / dx - 0.5f - j;
        float val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )          ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) + 
                     ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                     ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                     ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }
    
    /**
     * Interpolates the given velocityu
     * @param x
     * @param s
     * @return interpolated value
     */
    public float interpolateU( Tuple2f x, float[] s ) {
        int i,j;
        float wx, wy;
        i = (int) Math.floor( x.x / dx);
        j = (int) Math.floor( x.y / dx - 0.5 );
        wx = x.x / dx - i;
        wy = x.y / dx - 0.5f - j;
        float val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )          ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) + 
                     ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                     ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                     ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }
    
    /**
     * Interpolates the given velocityv
     * @param x
     * @param s
     * @return interpolated value
     */
    public float interpolateV( Tuple2f x, float[] s ) {
        int i,j;
        float wx, wy;
        i = (int) Math.floor( x.x / dx - 0.5 );
        j = (int) Math.floor( x.y / dx);
        wx = x.x / dx - 0.5f - i;
        wy = x.y / dx - j;
        float val = ( (i>=0 && j>=0 && i<=N+1 && j<=N+1 )          ? s[IX(i,j)]*(1-wx)*(1-wy) : 0 ) + 
                     ( (i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )     ? s[IX(i+1,j)]*wx*(1-wy) : 0 ) +
                     ( (i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )     ? s[IX(i,j+1)]*(1-wx)*wy : 0 ) +
                     ( (i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 ) ? s[IX(i+1,j+1)]*wx*wy : 0 );
        return val;
    }
    
    /** 
     * Performs a simple Forward Euler particle trace using the current velocity field 
     * @param x0
     * @param h
     * @param x1
     */
    public void traceParticle( Point2f x0, float h, Point2f x1 ) {
        traceParticle( x0, U0, h, x1 );        
    }
    
    /** 
     * Performs a simple particle trace using Rk3
     * x1 = x0 + h/6 *K1  + 2/3 * h * U(x0+0.5*h*k1) + h/6 * U(x0 - h * K1 + 2*h*K2)
     * K1=U(x0)
     * K2=U(x0+K1/2*h)
     * k3=U(x0-h*K1+2*K2*h)
     * @param x0
     * @param U
     * @param h
     * @param x1
     */
    private void traceParticle( Point2f x0, float[][] U, float h, Point2f x1 ) {
    	if (RK3.getValue()){
            Vector2f k1 = new Vector2f();
            Vector2f k2 = new Vector2f();
            Vector2f k3 = new Vector2f();
            Point2f x2 = new Point2f();
            Point2f x3 = new Point2f();
            x1.set( x0 );
            getVelocity(x1, U, k1);
            x2.x=x0.x + 0.5f*k1.x*h;
            x2.y=x0.y + 0.5f*k1.y*h;
            getVelocity(x2, U, k2);
            x3.x=x0.x - k1.x*h + 2*h*k2.x;
            x3.y=x0.y - k1.y*h + 2*h*k2.y;
            getVelocity(x3, U, k3);
            k1.scale(h/6.0f);
            k2.scale(2.0f*h/3.0f);
            k3.scale(h/6.0f);
            x1.add( k1 );
            x1.add( k2 );
            x1.add( k3 );
            //System.out.println("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL");
    	}else{//Forword eluer
    		  /** 
    	     * Performs a simple particle trace using Forward Euler
    	     * (Note that this could be something higher order or adative)
    	     * x1 = x0 + h * U(x0)
    	     * @param x0
    	     * @param U
    	     * @param h
    	     * @param x1
    	     */
    		Vector2f vel = new Vector2f();
            x1.set( x0 );
            getVelocity(x1, U, vel);
            vel.scale(h);            
            x1.add( vel );
    	}

    }
    
    private Point2f mouseX1 = new Point2f();
    
    private Point2f mouseX0 = new Point2f();
    
    private void addMouseForce( float[][] U, float dt ) {
                
        Vector2f f = new Vector2f();
        f.sub( mouseX1, mouseX0 );
        float d = mouseX0.distance(mouseX1);
        if ( d < 1e-6 ) return;
        f.scale( mouseForce.getFloatValue() );
        // go along the path of the mouse!
        Point2f x = new Point2f();
        int num = (int) (d/dx + 1);
        for ( int i = 0; i <= num; i++ ) {
            x.interpolate(mouseX0,mouseX1, (float)i / num );
            addForce( U, dt, x, f );
        }
        mouseX0.set( mouseX1 );
    }
    
    /**
     * Sets the mouse location in the fluid for doing dragging interaction
     * @param x0 previous location
     * @param x1 current location
     */
    public void setMouseMotionPos( Point2f x0, Point2f x1 ) {
        mouseX0.set( x0 );
        mouseX1.set( x1 );
    }
    
    /**
     * Adds a force at a given point in the provided velocity field.
     * @param U    velocity field
     * @param dt   time step
     * @param x    location
     * @param f    force
     */
    private void addForce( float[][] U, float dt, Tuple2f x, Tuple2f f ) {
    	addVelocityu( U[0], dt, x, f.x );
    	addVelocityv( U[1], dt, x, f.y );
    }
     
    /**
     * Adds the time step scaled amount to the provided scalar field at the specified location x in an interpolated manner.
     * @param S    scalar field
     * @param dt   time step
     * @param x    location
     * @param a    amount
     */
    private void addSource( float[] S, float dt, Tuple2f x, float a ) {
        int i = (int) Math.floor( x.x / dx - 0.5 );
        int j = (int) Math.floor( x.y / dx - 0.5 );
        float wx = x.x / dx - 0.5f - i;
        float wy = x.y / dx - 0.5f - j;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }
    
    /**
     * Adds the time step scaled amount to the provided horizontal velocity field at the specified location x in an interpolated manner.
     * @param S    scalar field
     * @param dt   time step
     * @param x    location
     * @param a    amount
     */
    private void addVelocityu( float[] S, float dt, Tuple2f x, float a ) {
    	int i = (int) Math.floor( x.x / dx );
        int j = (int) Math.floor( x.y / dx -0.5);
        float wx = x.x / dx - i ;
        float wy = x.y / dx - 0.5f - j;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }
    
    /**
     * Adds the time step scaled amount to the provided vertical velocity field at the specified location x in an interpolated manner.
     * @param S    scalar field
     * @param dt   time step
     * @param x    location
     * @param a    amount
     */
    private void addVelocityv( float[] S, float dt, Tuple2f x, float a ) {
        int i = (int) Math.floor( x.x / dx - 0.5 );
        int j = (int) Math.floor( x.y / dx );
        float wx = x.x / dx - 0.5f - i;
        float wy = x.y / dx - j ;
        if ( i>=0 && j>=0 && i<=N+1 && j<=N+1 )          S[IX(i,j)]     += (1-wx)*(1-wy)* dt * a;
        if ( i+1>=0 && j>=0 && i+1<=N+1 && j<=N+1 )      S[IX(i+1,j)]   += wx*(1-wy) * dt * a;
        if ( i>=0 && j+1>=0 && i<=N+1 && j+1<=N+1 )      S[IX(i,j+1)]   += (1-wx)*wy * dt * a;
        if ( i+1>=0 && j+1>=0 && i+1<=N+1 && j+1<=N+1 )  S[IX(i+1,j+1)] += wx*wy * dt * a;
    }
    /**
     * Gets the average temperature of the continuum.
     * (use this in computing bouyancy forces)
     * @return
     */
    public double getReferenceTemperature() {
        int count = 0;
        double referenceTemperature = 0;
        for ( int i = 1; i <= N; i++ ) {
            for ( int j = 1; j <= N; j++ ) {
                referenceTemperature += temperature0[IX(i,j)];
                count++;
            }
        }
        referenceTemperature /= count;
        return referenceTemperature;
    }
    
     
     
     
    // TODO: You may want to implement a bunch of methods here to diffuse, transport, fix boundaries, project, and step vector and scalar fields
    
    
    private void swap (float[] x, float[] x0){
    	int i, j;
		   for ( i=1 ; i<=N ; i++ ) {
			   for ( j=1 ; j<=N ; j++ ) {
				   float tmp = x0[IX(i,j)];
				   x0[IX(i,j)]= x[IX(i,j)];
				   x[IX(i,j)]=tmp;
			   }
		   }
    }
    

    private void diffuse(int b, float[] d, float[] d0, float diff, float dt){
        int i, j, k;
        int iter = iterations.getValue();
        float a = dt * diff * N * N;
                      
        for( k = 0; k < iter; k ++){
            
            for(i = 1 ; i <= N; i++){
                for( j = 1; j <= N; j++){
                    d[IX(i,j)] = (d0[IX(i,j)] + a*(d[IX(i-1, j)] + d[IX(i+1, j)] +d[IX(i, j - 1)] + d[IX(i,j + 1)])) / (1 + 4 * a);
                }
            }
            
        }
        set_bnd(N, b, d);   
    }
    
    private void diffuseU(int b, float[] d, float[] d0, float diff, float dt){
        int i, j, k;
        int iter = iterations.getValue();
        float a = dt * diff * N * N;
                      
        for( k = 0; k < iter; k ++){
            
            for(i = 1 ; i <= N; i++){
                for( j = 1; j <= N; j++){
                    d[IX(i,j)] = (d0[IX(i,j)] + a*(d[IX(i-1, j)] + d[IX(i+1, j)] +d[IX(i, j - 1)] + d[IX(i,j + 1)])) / (1 + 4 * a);
                }
            }
            
        }
        set_bndVelocityU(N, d);   
    }
    
    private void diffuseV(int b, float[] d, float[] d0, float diff, float dt){
        int i, j, k;
        int iter = iterations.getValue();
        float a = dt * diff * N * N;
                      
        for( k = 0; k < iter; k ++){
            
            for(i = 1 ; i <= N; i++){
                for( j = 1; j <= N; j++){
                    d[IX(i,j)] = (d0[IX(i,j)] + a*(d[IX(i-1, j)] + d[IX(i+1, j)] +d[IX(i, j - 1)] + d[IX(i,j + 1)])) / (1 + 4 * a);
                }
            }
            
        }
        set_bndVelocityV(N, d);  
    }
    
    
    
    private void advect( int b, float[] d, float[] d0, float[][] U, float dt){
        int i, j;
        float x,y;
         
        for (i = 1; i <= N; i++){
            for(j = 1; j <= N ; j++){
                
                x = (i + 0.5f) * dx;
                y = (j + 0.5f) * dx;
                Point2f x0 = new Point2f(x,y);
                Point2f x1 = new Point2f();
                traceParticle(x0, U, -dt, x1);
                d[IX(i,j)] = interpolate(x1,d0 );
            
                
            }
        }
        set_bnd(N,b,d);
    }
    
    private void advectU( int b, float[] d, float[] d0, float[][] U, float dt){
        int i, j;
        float x,y;
         
        for (i = 1; i <= N; i++){
            for(j = 1; j <= N ; j++){
                
                x = (i) * dx;
                y = (j + 0.5f) * dx;
                Point2f x0 = new Point2f(x,y);
                Point2f x1 = new Point2f();
                traceParticle(x0, U,-dt, x1);
                d[IX(i,j)] = interpolateU(x1,d0 );
            
                
            }
        }
        set_bndVelocityU(N,d);
    }
    
    private void advectV( int b, float[] d, float[] d0, float[][] U, float dt){
        int i, j;
        float x,y;
         
        for (i = 1; i <= N; i++){
            for(j = 1; j <= N ; j++){
                
                x = (i+0.5f) * dx;
                y = (j) * dx;
                Point2f x0 = new Point2f(x,y);
                Point2f x1 = new Point2f();
                traceParticle(x0, U,-dt, x1);
                d[IX(i,j)] = interpolateV(x1,d0 );
            
                
            }
        }
        set_bndVelocityV(N,d);
    }
    
    private void dens_step( float[][] U, float diff, float dt ){
        
        diffuse(0, temperature1, temperature0, diff, dt);
        swap(temperature1,temperature0);
        advect(0, temperature1, temperature0, U, dt);
        swap(temperature1,temperature0); 
    }
    

    private void set_bnd(int N, int b, float[] U){
        int i;
        for(i = 1; i <= N; i++){
            U[IX(0,i)] = b==1 ? -U[IX(1,i)] : U[IX(1,i)];
            U[IX(N+1,i)] = b==1 ? -U[IX(N,i)] : U[IX(N,i)];
            U[IX(i,0)] = b==2 ? -U[IX(i,1)] : U[IX(i,1)];
            U[IX(i,N+1)] = b==2 ? -U[IX(i,N)] : U[IX(i,N)];
        }
        
        U[IX(0,0)] = (float) (0.5 * (U[IX(1,0)] + U[IX(0,1)]));
        U[IX(0,N+1)] = (float) (0.5 * (U[IX(1,N+1)] + U[IX(0,N)]));
        U[IX(N+1,0)] = (float) (0.5 * (U[IX(N,0)] + U[IX(N+1,1)]));
        U[IX(N+1,N+1)] = (float) (0.5 * (U[IX(N,N+1)] + U[IX(N+1,N)]));
        
    }
    
    private void set_bndVelocityU(int N,float[] U){
        int j;
        for(j = 1; j <= N; j++){
            U[IX(1,j)] = 0.f;
            U[IX(N+1,j)] = 0.f;
        }
    }
    
    private void set_bndVelocityV(int N,float[] U){
        int i;
        for(i = 1; i <= N; i++){
            U[IX(i,1)] = 0.f;
            U[IX(i,N+1)] = 0.f;
        }
    }
    
    //INPUT MATRIX A diagonal is 4, -1 for the neighbour cells
    private void setA(){
    	for(int i =0;i<=N+1;i++){
    		for(int j =0;j<=N+1;j++){
    			diagonalA[IX(i,j)]=0;
                plusiA[IX(i,j)]=0;
                plusjA[IX(i,j)]=0;
    		}
    	}
    	for (int i = 1; i<= N; i++){
                    for(int j = 1; j <= N; j++){
                    diagonalA[IX(i,j)]=4f;
                    plusiA[IX(i,j)]=-1f;
                    plusjA[IX(i,j)]=-1f;
                    }
        }
    }

    //apply A multiply matrixA with a vector
    private void applyA(float[] x, float[] y)
    {
    	for(int j=0; j<=N+1; j++){ 
      	   for(int i=0; i<=N+1; i++){
      		 y[IX(i,j)]=0;
      	   }
      }
       for(int i=1; i<=N; i++){ 
    	   for(int j=1; j<=N; j++){
             y[IX(i,j)]=diagonalA[IX(i,j)]*x[IX(i,j)] + plusiA[IX(i-1,j)]*x[IX(i-1,j)]
                                          + plusiA[IX(i,j)]*x[IX(i+1,j)]
                                          + plusjA[IX(i,j-1)]*x[IX(i,j-1)]
                                          + plusjA[IX(i,j)]*x[IX(i,j+1)];
    	   }
       }
    }

    //Construct the MIC(0) preconditioner
       
    private void formPreconditioner()
    {
       float mic_parameter=0.97f;
       float e;
       for(int i=0; i<=N+1; i++){ 
      	   for(int j=0; j<=N+1; j++){
      		 preconditioner[IX(i,j)]=0f;
      	   }
      }
       
       for(int i=1;i<=N; i++) {
       		for(int j=1; j<=N; ++j){
             e=diagonalA[IX(i,j)] - ( plusiA[IX(i-1,j)]*preconditioner[IX(i-1,j)] ) * ( plusiA[IX(i-1,j)]*preconditioner[IX(i-1,j)] )
                              - ( plusjA[IX(i,j-1)]*preconditioner[IX(i,j-1)] ) *( plusjA[IX(i,j-1)]*preconditioner[IX(i,j-1)] )
                              - mic_parameter*( plusiA[IX(i-1,j)]*plusjA[IX(i-1,j)]*(preconditioner[IX(i-1,j)])*(preconditioner[IX(i-1,j)])
                                               +plusjA[IX(i,j-1)]*plusiA[IX(i,j-1)]*(preconditioner[IX(i,j-1)])*(preconditioner[IX(i,j-1)]) );
             preconditioner[IX(i,j)]=(float) (1f/Math.sqrt(e+1e-17f));
          }
       }
    }
    
    //apply preconditioner MIC(0)
    private void apply_preconditioner(float[] r, float[] z)
    {
       int i, j;
       float t,f;
       float[] q = new float[(N+2)*(N+2)];
       for(int k=0; k<N+1; k++){
      	   for(int l=0; l<N+1; l++){
      		 q[IX(k,l)]=0;
      	   }
      }
       // solve L*q=r
       for(i=1; i<=N; i++) {
    	   for(j=1; j<=N; j++){
             t=r[IX(i,j)] - plusiA[IX(i-1,j)]*preconditioner[IX(i-1,j)]*q[IX(i-1,j)]
                      - plusjA[IX(i,j-1)]*preconditioner[IX(i,j-1)]*q[IX(i,j-1)];
             q[IX(i,j)]=t*preconditioner[IX(i,j)];
    	   }
       }
       // solve L'*z=q
       //float[] z = new float[(N+2)*(N+2)];
       for(int k=0; k<N+1; k++){ 
      	   for(int l=0; l<N+1; l++){
      		 z[IX(k,l)]=0;
      	   }
      }
       for( i=N; i>=1; i--) {
    	   for( j=N; j>=1; j--){
             f=q[IX(i,j)] - plusiA[IX(i,j)]*preconditioner[IX(i,j)]*z[IX(i+1,j)]
                      - plusjA[IX(i,j)]*preconditioner[IX(i,j)]*z[IX(i,j+1)];
             z[IX(i,j)]=f*preconditioner[IX(i,j)];
    	   }
       }
    }
    
    //find the infinity norm of a vector
    private float infinityNorm(float[] h){
    	float infnorm =0;
    	for (int i = 1; i <= N; i++){
            for(int j = 1; j <= N; j++){
            	if(infnorm<=Math.abs(h[IX(i,j)])){
            		infnorm=Math.abs(h[IX(i,j)]);
            	}
   			}
   	}
		return infnorm;
    }
    
    //Caculate the dot product of two vector
    private float dotProduct(float[] h, float[] o){
    	float dotproduct=0;
    	for (int i = 1; i <= N; i++){
            for(int j = 1; j <= N; j++){
            	dotproduct += h[IX(i,j)]*o[IX(i,j)];
   			}
    	}
    	return dotproduct;
    }
    
    
    //pcg solver with Modified Incomplete Cholesky
    private void solvePressure(int maxits, float tol)
    {
       int its;
       float infnorm;
       for (int i = 1; i <= N; i++){
                for(int j = 1; j <= N; j++){
       				r[IX(i,j)]=div[IX(i,j)];
       				p[IX(i,j)]=0;
       			}
       	}
       infnorm = infinityNorm(r);

       //tol=tol*infinityNorm(r);
       
       if(infnorm==0){
         return;
       }
//
       apply_preconditioner(r,z);
       for (int i = 1; i <= N; i++){
           for(int j = 1; j <= N; j++){
  				s[IX(i,j)]=z[IX(i,j)];
  			}
  	}
       float theta=dotProduct(z,r);
       if(theta==0) {
    	   return;
       }
    	   
       for(its=0; its<100; its++){
          applyA(s, z);
          //System.out.println(infinityNorm(z));
          float alpha=theta/dotProduct(z,s);
          for (int i = 1; i <= N; i++){
              for(int j = 1; j <= N; j++){
     				p[IX(i,j)]+=alpha*s[IX(i,j)];
     			}
          }
          for (int i = 1; i <= N; i++){
              for(int j = 1; j <= N; j++){
     				r[IX(i,j)]-=alpha*z[IX(i,j)];
     			}
          }
          
          if(infinityNorm(r)<tol){
        	//System.out.printf("pressure converged to %g in %d iterations\n", infinityNorm(r), its);
        	 break;
          }

          apply_preconditioner(r,z);
          float thetanew=dotProduct(z,r);
          float beta=thetanew/theta;
          for (int i = 1; i <= N; i++){
              for(int j = 1; j <= N; j++){
     				s[IX(i,j)] = z[IX(i,j)] + beta*s[IX(i,j)];
     			}
          }
          theta=thetanew;

       }
       if (infinityNorm(r)>=tol && its==100){
           System.out.printf("Didn't converge in pressure solve (its=%d, tol=%g, |r|=%g)\n", its, tol, infinityNorm(r)); 
       }

       return;
    }
    private void project(int N, float[][] U){
        int i, j, k;
        float h;
        int iter = iterations.getValue();
        h = 1.0f / N;
        for (i = 1; i <= N; i++){
            for(j = 1; j <= N; j++){
                div[IX(i,j)] =  -h * (U[0][IX(i+1,j)] - U[0][IX(i,j)] + U[1][IX(i,j+1)] - U[1][IX(i,j)]);
            }
        }
        set_bnd(N,0,div);
        set_bnd(N, 0, p);
        
        // pcg solver 
        if (pcgsolver.getValue()){
            setA();
            formPreconditioner();
            solvePressure(100, 1e-20f);
            set_bnd(N, 0, p);
            float [] residual = new float[(N+2)*(N+2)];
            for(i=1; i<=N; i++){ 
         	   for(j=1; j<=N; j++){
                  residual[IX(i,j)]=div[IX(i,j)]-(diagonalA[IX(i,j)]*p[IX(i,j)] + plusiA[IX(i-1,j)]*p[IX(i-1,j)]
                                               + plusiA[IX(i,j)]*p[IX(i+1,j)]
                                               + plusjA[IX(i,j-1)]*p[IX(i,j-1)]
                                               + plusjA[IX(i,j)]*p[IX(i,j+1)]);
         	   }
            }
            System.out.println(infinityNorm(residual));
        }else{
            for(k = 0; k < iter; k++){
                for (i = 1; i<= N; i++){
                    for(j = 1; j <= N; j++){
                        p[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1, j)] + p[IX(i+1, j)] + p[IX(i,j-1)] + p[IX(i,j+1)])/4;
                    }
                }   
                set_bnd(N, 0, p); 
            }
            float [] residual = new float[(N+2)*(N+2)];
            for(i=1; i<=N; i++){ 
         	   for(j=1; j<=N; j++){
                  residual[IX(i,j)]=div[IX(i,j)]-(diagonalA[IX(i,j)]*p[IX(i,j)] + plusiA[IX(i-1,j)]*p[IX(i-1,j)]
                                               + plusiA[IX(i,j)]*p[IX(i+1,j)]
                                               + plusjA[IX(i,j-1)]*p[IX(i,j-1)]
                                               + plusjA[IX(i,j)]*p[IX(i,j+1)]);
         	   }
            }
            System.out.println(infinityNorm(residual));

        }//gs slover
            

        //velocity updates
        for(i = 1; i <= N; i++){
            for(j = 1; j <= N; j ++){
                
                U[0][IX(i,j)] -= (p[IX(i,j)] - p[IX(i-1,j)])/h;
                U[1][IX(i,j)] -= (p[IX(i,j)] - p[IX(i,j-1)])/h;
                //System.out.println(U[0][IX(i,j)]); 
            }
        }
        set_bndVelocityU(N, U[0]);
        set_bndVelocityV(N, U[1]);
    }
    
    private void vel_step( float visc, float dt){

    	//this.SwapU();
        diffuseU(1, U1[0], U0[0], visc, dt);
        diffuseV(2, U1[1], U0[1], visc, dt);

        project(N, U1);

        swap(U1[0],U0[0]);
        swap(U1[1],U0[1]);
        
        advectU(1, U1[0], U0[0], U0, dt);
        advectV(2, U1[1], U0[1], U0, dt);
        

        project(N, U1);
        swap(U1[0],U0[0]);
        swap(U1[1],U0[1]);
        
        
    }
    
    private float maxvelocity(float[] u, float[] v) {
		float max = 0;
    	for(int i=0;i<N+1;i++){
    		float step=(float) Math.sqrt(Math.abs(u[i]*u[i])+Math.abs(v[i]*v[i]));
    		if (max<=step){
    			max = step;
    		}
    	}
    	return max;
    }
    /**
     * Advances the state of the fluid by one time step
     */
    public void step() {
        
    	//calulate dt 
        //use cfl condition dh=1.0f/N
        //TO DO THIS WE HAVE TO FIND THE MAXq VELOCITY
        float umax=maxvelocity(U0[0],U0[1]);
        umax=(float) ((umax + Math.sqrt((1.0/N)*9.8)));
        dt=(float) (kcfl.getFloatValue()*(1.0/N)/umax);

        //dt=stepSize.getFloatValue();
        addMouseForce( U0, dt );        
        // TODO: use sources to change temperatures in the temperature0 scalar field

        for (Source s : sources){
            addSource(temperature0, dt, s.location, s.amount);
        }     
        // TODO: use temperature scalar field to apply buoyancy forces to the velocity field

        double T0 = getReferenceTemperature();
        float beta = buoyancy.getFloatValue();
        for(int i = 1; i <= N; i++){
        	for(int j = 1; j <= N; j++){
        		float Tk = temperature0[IX(i,j)] ;
        		U0[1][IX(i,j)] += beta*(T0-Tk)*9.8*dt;
        		//U0[1][IX(i,j)] += 9.8*dt/N;
        	}
        }
        
        
        // TODO: do the velocity step

       vel_step( viscosity.getFloatValue(), dt);
     

        // do the scalar step; 
       	// velocity step has swapped updated U1 to U0... 
        
        dens_step( U0, diffusion.getFloatValue(), dt);
        
        
        elapsed += dt;
    }

	private DoubleParameter kcfl = new DoubleParameter( "kCFL", 1, 1, 5 );    
    private DoubleParameter viscosity = new DoubleParameter( "viscosity", 1e-6, 1e-8, 1 );
    private DoubleParameter diffusion = new DoubleParameter( "diffusion", 1e-6, 1e-8, 1 );
    private DoubleParameter buoyancy = new DoubleParameter( "buoyancy", 0.1, -1, 1 );
    private IntParameter iterations = new IntParameter( "Gauss Seidel iterations", 30, 20, 200 );    
    private DoubleParameter mouseForce = new DoubleParameter( "mouse force", 1e2, 1, 1e3 );
    private BooleanParameter pcgsolver = new BooleanParameter("PCG solver", false);
    private BooleanParameter RK3 = new BooleanParameter("rk3", false);
    /** step size of the simulation */
    public DoubleParameter stepSize = new DoubleParameter( "step size", 0.1, 0.001, 1 );
    
    /**
     * Get the parameters for the fluid 
     * @return a control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder(new TitledBorder("Fluid properties"));
        vfp1.add( viscosity.getSliderControls(true) );
        vfp1.add( diffusion.getSliderControls(true) );
        vfp1.add( buoyancy.getSliderControls(false) );
        vfp1.add( mouseForce.getSliderControls(true) );
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder(new TitledBorder("Fluid solver properties"));
        vfp2.add( pcgsolver.getControls() );
        vfp2.add( RK3.getControls() );
        vfp2.add( kcfl.getSliderControls(true ) ); 
        vfp2.add( stepSize.getSliderControls(true ) ); 
        vfp2.add( iterations.getSliderControls() );
        vfp2.add( Nval.getSliderControls() );
        vfp.add( vfp1.getPanel() );
        vfp.add( vfp2.getPanel() );
        return vfp.getPanel();
    }
}