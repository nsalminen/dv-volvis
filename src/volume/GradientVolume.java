/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel and modified by Anna Vilanova 
 * 
 * 
 */

//////////////////////////////////////////////////////////////////////
///////////////// CONTAINS FUNCTIONS TO BE IMPLEMENTED ///////////////
//////////////////////////////////////////////////////////////////////
public class GradientVolume {

	
    
//Do NOT modify this attributes
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Float[] laplace;
    
    Volume volume;
    double maxmag;
    
//If needed add new attributes here:
    private float maxLaplace;


    //Do NOT modify this function
    // 
    // Computes the gradient of the volume attribute and save it into the data attribute
    // This is a lengthy computation and is performed only once (have a look at the constructor GradientVolume) 
    //
    private void compute() {

        for (int i=0; i<data.length; i++) {
            data[i] = zero;
            laplace[i] = 0.0f;
        }
       
        for (int z=1; z<dimZ-1; z++) {
            for (int y=1; y<dimY-1; y++) {
                for (int x=1; x<dimX-1; x++) {
                    float gx = (volume.getVoxel(x+1, y, z) - volume.getVoxel(x-1, y, z))/2.0f;
                    float gy = (volume.getVoxel(x, y+1, z) - volume.getVoxel(x, y-1, z))/2.0f;
                    float gz = (volume.getVoxel(x, y, z+1) - volume.getVoxel(x, y, z-1))/2.0f;
                    VoxelGradient grad = new VoxelGradient(gx, gy, gz);
                    setGradient(x, y, z, grad);
                }
            }
        }
        
        for (int z=2; z<dimZ-2; z++) {
            for (int y=2; y<dimY-2; y++) {
                for (int x=2; x<dimX-2; x++) {         
                    float l = volume.getVoxel(x-1, y, z) + volume.getVoxel(x+1, y, z) 
                            +  volume.getVoxel(x, y + 1, z) + volume.getVoxel(x, y - 1, z)
                             +  volume.getVoxel(x, y, z + 1) + volume.getVoxel(x, y, z - 1) 
                            - 6 * volume.getVoxel(x, y, z); 
                   
                    setLaplace(x, y, z, l);
                }
            }
        }
        maxmag=calculateMaxGradientMagnitude();
        maxLaplace = calculateMaxLaplace();
     }
    	
    
//////////////////////////////////////////////////////////////////////
///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
//////////////////////////////////////////////////////////////////////
//This function linearly interpolates gradient vector g0 and g1 given the factor (t) 
//the resut is given at result. You can use it to tri-linearly interpolate the gradient 
    
    // This function linearly interpolates gradient vector g0 and g1 given the factor (t), and returns
    // the value through the result variable
    private void interpolate(VoxelGradient g0, VoxelGradient g1, float factor, VoxelGradient result) {        
        result.x = factor * g0.x + (1 - factor) * g1.x;
        result.y = factor * g0.y + (1 - factor) * g1.y;
        result.z = factor * g0.z + (1 - factor) * g1.z;
        result.mag = (float) Math.sqrt(result.x * result.x + result.y * result.y + result.z * result.z);
    }
    
    // This function linearly interpolates gradient vector g0 and g1 given the factor (t), and returns
    // the value in a new VoxelGradient object
    private VoxelGradient interpolate(VoxelGradient g0, VoxelGradient g1, float factor) {     
        VoxelGradient result = new VoxelGradient();
        result.x = factor * g0.x + (1 - factor) * g1.x;
        result.y = factor * g0.y + (1 - factor) * g1.y;
        result.z = factor * g0.z + (1 - factor) * g1.z;
        result.mag = (float) Math.sqrt(result.x * result.x + result.y * result.y + result.z * result.z);
        return result;
    }
	
//////////////////////////////////////////////////////////////////////
///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
//////////////////////////////////////////////////////////////////////    
    // Returns the linearly interpolated gradient for position coord[]
    public VoxelGradient getGradient(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }
        
        // Compute the rounded up/down coordinate values
        double xF = Math.floor(coord[0]);
        double xC = Math.ceil(coord[0]);
        double yF = Math.floor(coord[1]);
        double yC = Math.ceil(coord[1]);
        double zF = Math.floor(coord[2]);
        double zC = Math.ceil(coord[2]);
        
        float dx = (float)(coord[0] - xF);
        float dy = (float)(coord[1] - yF);
        float dz = (float)(coord[2] - zF);
        
        // Interpolate along the x-axis
        VoxelGradient c00 = interpolate(getGradient((int) xF, (int) yF, (int) zF),
                                        getGradient((int) xC, (int) yF, (int) zF), 1.f - dx);
        VoxelGradient c01 = interpolate(getGradient((int) xF, (int) yF, (int) zC),
                                        getGradient((int) xC, (int) yF, (int) zC), 1.f - dx);
        VoxelGradient c10 = interpolate(getGradient((int) xF, (int) yC, (int) zF),
                                        getGradient((int) xC, (int) yC, (int) zF), 1.f - dx);
        VoxelGradient c11 = interpolate(getGradient((int) xF, (int) yC, (int) zC),
                                        getGradient((int) xC, (int) yC, (int) zC), 1.f - dx);
        
        // Interpolate along the y-axis
        VoxelGradient c0 = interpolate(c00,c10, 1.f - dy);
        VoxelGradient c1 = interpolate(c01,c11, 1.f - dy);
        
        // Interpolate along the z-axis
        VoxelGradient c = interpolate(c0, c1, 1.f - dz);
        
        return c;
    }
    
    
     public float getLaplace(double[] coord) {
        if (coord[0] < 1 || coord[0] > (dimX-3) || coord[1] < 1 || coord[1] > (dimY-3)
                || coord[2] < 1 || coord[2] > (dimZ-3)) {
            return 0;
        }
        
        // Compute the rounded up/down coordinate values
        double xF = Math.floor(coord[0]);
        double xC = Math.ceil(coord[0]);
        double yF = Math.floor(coord[1]);
        double yC = Math.ceil(coord[1]);
        double zF = Math.floor(coord[2]);
        double zC = Math.ceil(coord[2]);
        
        float dx = (float)(coord[0] - xF);
        float dy = (float)(coord[1] - yF);
        float dz = (float)(coord[2] - zF);
        
        // Interpolate along the x-axis
        float c00 = getLaplace((int) xF, (int) yF, (int) zF) * (1-dx) + 
                    getLaplace((int) xC, (int) yF, (int) zF) * dx;
        float c01 = getLaplace((int) xF, (int) yF, (int) zC) * (1-dx) + 
                    getLaplace((int) xC, (int) yF, (int) zC) * dx;
        float c10 = getLaplace((int) xF, (int) yC, (int) zF) * (1-dx) + 
                    getLaplace((int) xC, (int) yC, (int) zF) * dx;
        float c11 = getLaplace((int) xF, (int) yC, (int) zC) * (1-dx) + 
                    getLaplace((int) xC, (int) yC, (int) zC) * dx;
                                
   
        // Interpolate along the y-axis
        float c0 = c00 * (1 - dy) + c10 * dy;
        float c1 = c01* (1 - dy) + c11 * dy;
        
        // Interpolate along the z-axis
        float c = c0 * ( 1.f - dz) +  c1 * dz;
        return c;
    }
    
    
    //Do NOT modify this function
    public VoxelGradient getGradientNN(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        return getGradient(x, y, z);
    }
    
    //Returns the maximum gradient magnitude
    //The data array contains all the gradients, in this function you have to return the maximum magnitude of the vectors in data[] 
    //Do NOT modify this function
    private double calculateMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private float calculateMaxLaplace() {
        if (maxLaplace >= 0.f) {
            return maxLaplace;
        } else {
            float max = laplace[0];
            for (int i=0; i<laplace.length; i++) {
                max = laplace[i] > max ? laplace[i] : max;
            }   
            maxLaplace = max;
            return maxLaplace;
        }
    }
        
    //Do NOT modify this function
    public double getMaxGradientMagnitude()
    {
        return this.maxmag;
    }
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
	
	//Do NOT modify this function
	public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        laplace = new Float[dimX * dimY * dimZ];
        maxmag = -1.0;
        compute();
    }

	//Do NOT modify this function
	public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }
        
    public float getLaplace(int x, int y, int z) {
        return laplace[x + dimX * (y + dimY * z)];
    }

  
  
    //Do NOT modify this function
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }
    
    public void setLaplace(int x, int y, int z, float value) {
        laplace[x + dimX * (y + dimY * z)] = value;
    }

    //Do NOT modify this function
    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }
    
    //Do NOT modify this function
    public VoxelGradient getVoxel(int i) {
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

    public int getMaxLaplace() {
        return (int) maxLaplace;
    }

}
