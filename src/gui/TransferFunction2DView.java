/*
 * Anna Vilanova: Basic user interface for the raycaster. NO MODIFICATION NEEDED 
 *
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import volvis.TransferFunction2D;

/**
 *
 * @author michel and modified by Anna Vilanova 
 * 
 * NO IMPLEMENTATION NEEDED FOR THE BASIC ASSIGNMENTME 
 * YOu need to modify it if you want to change the interface
 */

public class TransferFunction2DView extends javax.swing.JPanel {

    TransferFunction2DEditor ed;
    private final int DOTSIZE = 8;
    public Ellipse2D.Double baseControlPoint, radiusControlPoint;
    public Ellipse2D.Double minMagnitudeControlPoint, maxMagnitudeControlPoint;
    boolean selectedBaseControlPoint, selectedRadiusControlPoint;
    boolean selectedMinMagnitudeControlPoint, selectedMaxMagnitudeControlPoint;
    
    
    /**
     * Creates new form TransferFunction2DView
     * @param ed
     */
    public TransferFunction2DView(TransferFunction2DEditor ed) {
        initComponents();
        
        this.ed = ed;
        selectedBaseControlPoint = false;
        selectedRadiusControlPoint = false;
        selectedMinMagnitudeControlPoint = false;
        selectedMaxMagnitudeControlPoint = false;
        addMouseMotionListener(new TriangleWidgetHandler());
        addMouseListener(new SelectionHandler());
    }
    
    @Override
    public void paintComponent(Graphics g) {

        Graphics2D g2 = (Graphics2D) g;

        int w = this.getWidth();
        int h = this.getHeight();
        g2.setColor(Color.white);
        g2.fillRect(0, 0, w, h);
        
        double maxHistoMagnitude = ed.histogram[0];
        for (int i = 0; i < ed.histogram.length; i++) {
            maxHistoMagnitude = ed.histogram[i] > maxHistoMagnitude ? ed.histogram[i] : maxHistoMagnitude;
        }
        
        double binWidth = (double) w / (double) ed.xbins;
        double binHeight = (double) h / (double) ed.ybins;
        maxHistoMagnitude = Math.log(maxHistoMagnitude);
        
        
        for (int y = 0; y < ed.ybins; y++) {
            for (int x = 0; x < ed.xbins; x++) {
                if (ed.histogram[y * ed.xbins + x] > 0) {
                    int intensity = (int) Math.floor(255 * (1.0 - Math.log(ed.histogram[y * ed.xbins + x]) / maxHistoMagnitude));
                    g2.setColor(new Color(intensity, intensity, intensity));
                    g2.fill(new Rectangle2D.Double(x * binWidth, h - (y * binHeight), binWidth, binHeight));
                }
            }
        }
        
        int ypos = h;
        int xpos = (int) (ed.tf2D.baseIntensity * binWidth);
        g2.setColor(Color.black);
        baseControlPoint = new Ellipse2D.Double(xpos - DOTSIZE / 2, ypos - DOTSIZE, DOTSIZE, DOTSIZE);
        g2.fill(baseControlPoint);
        g2.drawLine(xpos, ypos, xpos - (int) (ed.tf2D.radius * binWidth ), 0);
        g2.drawLine(xpos, ypos, xpos + (int) (ed.tf2D.radius * binWidth ), 0);
        radiusControlPoint = new Ellipse2D.Double(xpos + (ed.tf2D.radius * binWidth ) - DOTSIZE / 2,  0, DOTSIZE, DOTSIZE);
        g2.fill(radiusControlPoint);
        
        int maxY = (int) (h - ed.tf2D.maxMagnitude / ed.maxGradientMagnitude  * h);
        int minY = (int) (h - ed.tf2D.minMagnitude  / ed.maxGradientMagnitude  * h);
        maxY = maxY > h ? h : h <= 0 ? 1 :  maxY;
        minY = minY >= h ? h - 1 : h < 0 ? 0 : minY;
        g2.setColor(Color.red);
        maxMagnitudeControlPoint = new Ellipse2D.Double(0, maxY, DOTSIZE, DOTSIZE);
        g2.fill(maxMagnitudeControlPoint);
        g2.drawLine(0, maxY, w, maxY);
        
        g2.setColor(Color.blue);
        minMagnitudeControlPoint = new Ellipse2D.Double(0, minY - DOTSIZE, DOTSIZE, DOTSIZE);
        g2.fill(minMagnitudeControlPoint);
        g2.drawLine(0, minY, w, minY);
    }
    
    
    private class TriangleWidgetHandler extends MouseMotionAdapter {

        @Override
        public void mouseMoved(MouseEvent e) {
            
            if (baseControlPoint.contains(e.getPoint()) || radiusControlPoint.contains(e.getPoint())
                    || minMagnitudeControlPoint.contains(e.getPoint()) || maxMagnitudeControlPoint.contains(e.getPoint())) {
                setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            } else {
                setCursor(Cursor.getDefaultCursor());
            }
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            if (selectedBaseControlPoint || selectedRadiusControlPoint || selectedMinMagnitudeControlPoint || selectedMaxMagnitudeControlPoint) {
                Point dragEnd = e.getPoint();
                
                if (selectedBaseControlPoint) {
                    // restrain to horizontal movement
                    dragEnd.setLocation(dragEnd.x, baseControlPoint.getCenterY());
                } else if (selectedRadiusControlPoint) {
                    // restrain to horizontal movement and avoid radius getting 0
                    dragEnd.setLocation(dragEnd.x, radiusControlPoint.getCenterY());
                    if (dragEnd.x - baseControlPoint.getCenterX() <= 0) {
                        dragEnd.x = (int) (baseControlPoint.getCenterX() + 1);
                    }
                } else if (selectedMinMagnitudeControlPoint) {
                    dragEnd.setLocation(minMagnitudeControlPoint.getCenterX(), dragEnd.y);
                    if (dragEnd.y <= maxMagnitudeControlPoint.getCenterY() + DOTSIZE) {
                        dragEnd.y  = (int) minMagnitudeControlPoint.getCenterY() + DOTSIZE;
                    }
                } else if (selectedMaxMagnitudeControlPoint) {
                    dragEnd.setLocation(maxMagnitudeControlPoint.getCenterX(), dragEnd.y);
                    if (dragEnd.y >= minMagnitudeControlPoint.getCenterY() - DOTSIZE) {
                        dragEnd.y  = (int) maxMagnitudeControlPoint.getCenterY() - DOTSIZE;
                    }
                }
                if (dragEnd.x < 0) {
                    dragEnd.x = 0;
                }
                if (dragEnd.x >= getWidth()) {
                    dragEnd.x = getWidth() - 1;
                }
                if (dragEnd.y < 0) {
                    dragEnd.y = 0;
                }
                if (dragEnd.y >= getWidth()) {
                    dragEnd.y = getWidth() - 1;
                }
                double w = getWidth();
                double h = getHeight();
                double binWidth = (double) w / (double) ed.xbins;
                if (selectedBaseControlPoint) {
                    ed.tf2D.baseIntensity = (short) (dragEnd.x / binWidth);
                } else if (selectedRadiusControlPoint) {
                    ed.tf2D.radius = (dragEnd.x - (ed.tf2D.baseIntensity * binWidth))/binWidth;
                } else if (selectedMinMagnitudeControlPoint) {
                    ed.tf2D.minMagnitude = (h - dragEnd.y) / h * ed.maxGradientMagnitude;
                } else if (selectedMaxMagnitudeControlPoint) {
                    ed.tf2D.maxMagnitude = (h - dragEnd.y) / h * ed.maxGradientMagnitude;
                }
                ed.setSelectedInfo();
                
                repaint();
            } 
        }

    }
    
    
    private class SelectionHandler extends MouseAdapter {
        @Override
        public void mousePressed(MouseEvent e) {
            if (baseControlPoint.contains(e.getPoint())) {
                selectedBaseControlPoint = true;
            } else if (radiusControlPoint.contains(e.getPoint())) {
                selectedRadiusControlPoint = true;
            } else if (minMagnitudeControlPoint.contains(e.getPoint())) {
                selectedMinMagnitudeControlPoint = true;
            } else if (maxMagnitudeControlPoint.contains(e.getPoint())) {
                selectedMaxMagnitudeControlPoint = true;
            } else {
                selectedRadiusControlPoint = false;
                selectedBaseControlPoint = false;
                selectedMinMagnitudeControlPoint = false;
                selectedMaxMagnitudeControlPoint = false;
            }
        }
        
        @Override
        public void mouseReleased(MouseEvent e) {
            selectedRadiusControlPoint = false;
            selectedBaseControlPoint = false;
            selectedMinMagnitudeControlPoint = false;
            selectedMaxMagnitudeControlPoint = false;
            ed.changed();
            repaint();
        }
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
