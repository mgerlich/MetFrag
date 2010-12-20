/*
*
* Copyright (C) 2009-2010 IPB Halle, Sebastian Wolf
*
* Contact: swolf@ipb-halle.de
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/
package de.ipbhalle.metfrag.tools.renderer;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.RendererModel.ColorHash;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.font.IFontManager;
import org.openscience.cdk.renderer.generators.AtomNumberGenerator;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator.AtomRadius;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator.ColorByType;
import org.openscience.cdk.renderer.generators.BasicBondGenerator;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.IGenerator;
import org.openscience.cdk.renderer.generators.IGeneratorParameter;
import org.openscience.cdk.renderer.generators.RadicalGenerator;
import org.openscience.cdk.renderer.generators.RingGenerator;
import org.openscience.cdk.renderer.generators.BasicAtomGenerator.KekuleStructure;
import org.openscience.cdk.renderer.visitor.AWTDrawVisitor;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;


/**
* Example code for implementing a scrolling panel.
*
* @author maclean
*
*/
public class StructureRenderer extends JFrame {
    
    public class MoleculePanel extends JPanel {
        
        private int initialWidth;
        private int initialHeight;
        private AtomContainerRenderer renderer;
        private IAtomContainer atomContainer; 
        private boolean isNew;
        
        public MoleculePanel(IAtomContainer atomContainer) {
            
        	
        	this.atomContainer = atomContainer;
            
            this.initialWidth = 300;
            this.initialHeight = 300;
            
            this.setPreferredSize(
                    new Dimension(this.initialWidth + 10, this.initialHeight + 10));
            this.setBackground(Color.WHITE);
            this.setBorder(BorderFactory.createRaisedBevelBorder());
            
            List<IGenerator<IAtomContainer>> generators = new ArrayList<IGenerator<IAtomContainer>>();
            generators.add(new BasicSceneGenerator());
            generators.add(new BasicBondGenerator());
            generators.add(new BasicAtomGenerator());
            generators.add(new AtomNumberGenerator());
            generators.add(new RingGenerator());
            generators.add(new RadicalGenerator());
            
                                      
            
            IFontManager fm = new AWTFontManager();
            this.renderer = new AtomContainerRenderer(generators, fm); 
            RendererModel rm = renderer.getRenderer2DModel();
//            List<IGeneratorParameter<?>> parameterList = rm.getRenderingParameters();
//            for (IGeneratorParameter<?> parameter : parameterList) {
//				System.out.println(parameter.getClass().getName() + ": " +  parameter.getValue());
//			}
            
            
//            rm.set(ShowAromaticity.class, true);
            rm.set(KekuleStructure.class, true); 
            rm.set(AtomNumberGenerator.Offset.class, new javax.vecmath.Vector2d(15,0));
                        
            this.isNew = true;
        }
        
 
        public MoleculePanel(IAtomContainer atomContainer, List<Integer> atomsToHighlight) throws CDKException, IOException, CloneNotSupportedException {

        	this.atomContainer = atomContainer;           
            this.initialWidth = 300;
            this.initialHeight = 300;
            
            this.setPreferredSize(new Dimension(this.initialWidth + 10, this.initialHeight + 10));
            this.setBackground(Color.WHITE);
            this.setBorder(BorderFactory.createRaisedBevelBorder());
            
            List<IGenerator<IAtomContainer>> generators = new ArrayList<IGenerator<IAtomContainer>>();
            
            generators.add(new BasicSceneGenerator());
            generators.add(new BasicBondGenerator());
            generators.add(new BasicAtomGenerator());
//            generators.add(new RingGenerator());
             
            IFontManager fm = new AWTFontManager();
            this.renderer = new AtomContainerRenderer(generators, fm); 
            RendererModel rm = renderer.getRenderer2DModel();
   
            List<IBond> bondsToHighlight = new ArrayList<IBond>();
            Map<IChemObject, Color> colorMap = new HashMap<IChemObject, Color>();
            
            List<IAtom> atomsMatched = new ArrayList<IAtom>();
			for (Integer integer : atomsToHighlight) {
				atomsMatched.add(atomContainer.getAtom(integer));			   
			}
            
            for (IAtom atom : atomsMatched) {
            	for (IAtom atom2 : atomsMatched) {
            		
                	IBond bond = atomContainer.getBond(atom, atom2);
                	if(bond!= null)
                		bondsToHighlight.add(bond);
                }            	
            }
            for (IBond bond : bondsToHighlight) {
                colorMap.put(bond, new Color(0, 255, 0));
            }
            rm.getParameter(ColorHash.class).setValue(colorMap);
            
//            List<IGeneratorParameter<?>> parameterList = rm.getRenderingParameters();
//	        for (IGeneratorParameter<?> parameter : parameterList) {
//	        	System.out.println(parameter.getClass().getName() + ": " +  parameter.getValue());
//			}
                        
            this.isNew = true;
        }
        
        
        
        
        public void paint(Graphics g) {
            super.paint(g);
            
            
            Rectangle drawArea =
                new Rectangle(0, 0, this.initialWidth, this.initialHeight);
            
            if (this.isNew) {
                this.renderer.setup(atomContainer, drawArea);
                this.isNew = false;
            }
            
            Rectangle diagramRectangle =
                this.renderer.calculateDiagramBounds(atomContainer);
            
            Rectangle result = renderer.shift(drawArea, diagramRectangle);
            this.setPreferredSize(new Dimension(result.width, result.height));
            this.revalidate();
            
            this.renderer.paint(this.atomContainer, new AWTDrawVisitor((Graphics2D) g), drawArea, true);
        }
    }
    
    /**
     * Instantiates a new structure renderer. Just an example.
     */
    public StructureRenderer() {
        IMolecule chain = MoleculeFactory.makeAdenine();
                
        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule((IMolecule)chain);
        
        try {
            sdg.generateCoordinates();
            MoleculePanel molPanel = new MoleculePanel(sdg.getMolecule());
            this.add(new JScrollPane(molPanel));
        } catch (Exception e) {}
        
        this.pack();
        this.setVisible(true);
    }
    
    
    /**
     * Instantiates a new structure renderer.
     * Displays the structure in a window
     * 
     * @param original the original
     * @param name the name
     */
    public StructureRenderer(IAtomContainer original, String name) {
    	
    	IMolecule mol = new Molecule(original);
    	
    	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(mol);
        
        try {
            sdg.generateCoordinates();
            MoleculePanel molPanel = new MoleculePanel(sdg.getMolecule());
            this.add(new JScrollPane(molPanel));
        } catch (Exception e) {}
        
        this.pack();
        this.setVisible(true);
    }
    
    
    
    public StructureRenderer(IAtomContainer original, List<Integer> atomsTohighlight, String name) {
    	
    	IMolecule mol = new Molecule(original);
    	
    	StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(mol);
        
        try {
            sdg.generateCoordinates();
            MoleculePanel molPanel = new MoleculePanel(sdg.getMolecule(), atomsTohighlight);
            this.add(new JScrollPane(molPanel));
        } catch (Exception e) {
        	System.err.println("Error: " + e.getMessage() + "\n\n");
        	e.printStackTrace();
        }
        
        this.pack();
        this.setVisible(true);
    }


}
