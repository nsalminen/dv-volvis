/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel modified by Anna
 */

//////////////////////////////////////////////////////////////////////
///////////////// CONTAINS FUNCTIONS TO BE IMPLEMENTED ///////////////
//////////////////////////////////////////////////////////////////////

public class Volume {
    

    //Do NOT modify these attributes
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;

    // Do NOT modify this function
    // This function returns the nearest neighbour given a position in the volume given by coord.
    public short getVoxelNN(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-1) || coord[1] < 0 || coord[1] > (dimY-1)
                || coord[2] < 0 || coord[2] > (dimZ-1)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.round(coord[0]); 
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
    
        return getVoxel(x, y, z);
    }
        

    //Do NOT modify this function
    //This function linearly interpolates the value g0 and g1 given the factor (t) 
    //the result is returned. It is used for the tri-linearly interpolation the values 
    private float interpolate(float g0, float g1, float factor) {
        float result = (1 - factor)*g0 + factor*g1;
        return result; 
    }
             
    //Do NOT modify this function
    // This function returns the trilinear interpolated value of the position given by  position coord.
    public float getVoxelLinearInterpolate(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.floor(coord[0]); 
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        float fac_x = (float) coord[0] - x;
        float fac_y = (float) coord[1] - y;
        float fac_z = (float) coord[2] - z;

        float t0 = interpolate(getVoxel(x, y, z), getVoxel(x+1, y, z), fac_x);
        float t1 = interpolate(getVoxel(x, y+1, z), getVoxel(x+1, y+1, z), fac_x);
        float t2 = interpolate(getVoxel(x, y, z+1), getVoxel(x+1, y, z+1), fac_x);
        float t3 = interpolate(getVoxel(x, y+1, z+1), getVoxel(x+1, y+1, z+1), fac_x);
        float t4 = interpolate(t0, t1, fac_y);
        float t5 = interpolate(t2, t3, fac_y);
        float t6 = interpolate(t4, t5, fac_z);
        
        return t6; 
    }


    float a= 0.0f; // global variable that defines the value of a used in cubic interpolation.
    // you need to chose the right value
        
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
        
    // Function that computes the weights for one of the 4 samples involved in the 1D interpolation 
    private float weight (float x)//, Boolean one_two_sample) TODO no clue what the one_two_sample is
    {
        float a = -0.5f;
        x = (x < 0) ? -x : x;
        float result = 0;
        if (x < 1.0) {
            result = (a + 2) * x * x * x - (a + 3) * x * x + 1;
        } else if (x < 2.0) {
             result = a * x * x * x - 5 * a * x * x + 8 * a * x - 4 * a;
        } else {
            result = 0;
        }
       
        return (float)result; 
   }
    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // Function that computes the 1D cubic interpolation. g0,g1,g2,g3 contain the values of the voxels that we want to interpolate
    // factor contains the distance from the value g1 to the position we want to interpolate to.
    // We assume the out of bounce checks have been done earlier
    
    private float cubicinterpolate(float g0, float g1, float g2, float g3, float factor) {
       
        float dx0 = 1 + factor;
        float dx1 = factor;
        float dx2 = 1 - factor;
        float dx3 = 2 - factor;
        
        float result = weight(dx0) * g0 + weight(dx1) * g1 + weight(dx2) * g2 + weight(dx3) * g3;
                            
        return result; 
    }
        
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // 2D cubic interpolation implemented here. We do it for plane XY. Coord contains the position.
    // We assume the out of bounce checks have been done earlier
    private float bicubicinterpolateXY(double[] coord,int z) {
        
        float x = (float)coord[0];
        float y = (float)coord[1];
        //get coords for points around coord
        int x1 = (int) Math.floor(coord[0]); 
        int y1 = (int) Math.floor(coord[1]);        
        
        int x0 = x1 - 1;
        int y0 = y1 - 1;
     
        int x2 = x1 + 1;
        int y2 = y1 + 1;
        
        int x3 = x1 + 2;
        int y3 = y1 + 2;
        
        //TODO remove code duplication
        float t0 = cubicinterpolate(
                getVoxel(x0, y0, z),
                getVoxel(x1, y0, z),
                getVoxel(x2, y0, z),
                getVoxel(x3, y0, z),
                x - x1);      
        float t1 = cubicinterpolate(
                getVoxel(x0, y1, z),
                getVoxel(x1, y1, z),
                getVoxel(x2, y1, z),
                getVoxel(x3, y1, z),
                x - x1); 
        float t2 = cubicinterpolate(
                getVoxel(x0, y2, z),
                getVoxel(x1, y2, z),
                getVoxel(x2, y2, z),
                getVoxel(x3, y2, z),
                x - x1); 
        float t3 = cubicinterpolate(
                getVoxel(x0, y3, z),
                getVoxel(x1, y3, z),
                getVoxel(x2, y3, z),
                getVoxel(x3, y3, z),
                x - x1); 
        
        float result = cubicinterpolate(t0, t1, t2, t3, y - y1);
                            
        return result; 

    }
            
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // 3D cubic interpolation implemented here given a position in the volume given by coord.
    
    public float getVoxelTriCubicInterpolate(double[] coord) {
        if (coord[0] < 1 || coord[0] > (dimX-3) || coord[1] < 1 || coord[1] > (dimY-3)
                || coord[2] < 1 || coord[2] > (dimZ-3)) {
            return 0;
        }
        float z = (float)coord[2];
        
        int z1 = (int) Math.floor(coord[2]); 
        int z0 = z1 - 1;
        int z2 = z1 + 1;
        int z3 = z1 + 2;
        
        float t0 = bicubicinterpolateXY(coord, z0);
        float t1 = bicubicinterpolateXY(coord, z1);
        float t2 = bicubicinterpolateXY(coord, z2);
        float t3 = bicubicinterpolateXY(coord, z3);
                     
        float result = cubicinterpolate(t0, t1, t2, t3, z - z1);
        result = result < 0 ? 0 : result > 255 ? 255 : result;                
        return result; 
    }


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

	
    //Do NOT modify this function
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
    //Do NOT modify this function
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
    //Do NOT modify this function
    public short getVoxel(int x, int y, int z) {
    	int i = x + dimX*(y + dimY * z);
        return data[i];
    }
    
    //Do NOT modify this function
    public void setVoxel(int x, int y, int z, short value) {
    	int i = x + dimX*(y + dimY * z);
        data[i] = value;
    }
    
	//Do NOT modify this function
    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
    //Do NOT modify this function
    public short getVoxel(int i) {
        return data[i];
    }
    
	//Do NOT modify this function
    public int getDimX() {
        return dimX;
    }
    
    //Do NOT modify this function
    public int getDimY() {
        return dimY;
    }
    
    //Do NOT modify this function
    public int getDimZ() {
        return dimZ;
    }

    //Do NOT modify this function
    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }
    
    //Do NOT modify this function
    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
    //Do NOT modify this function
    public int[] getHistogram() {
        return histogram;
    }
    
    //Do NOT modify this function
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
}
