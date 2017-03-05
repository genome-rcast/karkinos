/*
Copyright Hiroki Ueda

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package jp.ac.utokyo.rcast.karkinos.graph;


//package com.pilotsw.pw.chart;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.Effect3D;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.ui.RectangleAnchor;
import org.jfree.chart.axis.AxisState;

import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Point2D;
import org.jfree.ui.Size2D;
import java.awt.geom.Line2D;
import java.awt.Shape;
import org.jfree.ui.RectangleEdge;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.axis.CategoryAnchor;
import org.jfree.chart.axis.CategoryLabelWidthType;
import org.jfree.chart.axis.CategoryTick;
import org.jfree.chart.axis.CategoryLabelPosition;
import org.jfree.chart.axis.Tick;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.entity.TickLabelEntity;
import org.jfree.data.category.CategoryDataset;
import org.jfree.text.TextBlock;
import org.jfree.text.TextLine;
import java.util.List;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
/*
 CategoryAxisSkipLabels skips category axis labels. Essentially this class extends
 the existing CategoryAxis class in JFreeChart, and overrides
 several of its methods (createLabel(), refreshTicks(), drawCategoryLabels()). 
 Occasionally I broke up the original code form JFreeChart::CategoryAxis() 
 into separate methods to make this fit together better and simplify testing. 

 It is important to realize that the original code, CategoryAxis, calculated the 
 positions of the label and their corresponding tick mark as they where drawn. 
 This technique makes it rather difficult to skip labels, certainly in a 
 well-balanced and efficient way. This class takes a different approach to the 
 problem of skipping labels. CategoryAxisSkipLabels calculates, and remembers, 
 given the current geometries, where each label will be draw. After all the 
 label positions are figured out and stored in the LabelsDefs list 
 CategoryAxisSkipLabels then decides which labels should be drawn.

 The central skipping algorithm is implemented in the removeOverlays() method. 
 This method accepts a plotting range, a start and end position of labels, and 
 attempts to draw the label at the mid point label of that range. If the label 
 does not collide with its marked for display neighbors it is marked for 
 display. In either case (marked or not) the current region is split and the 
 removeOverlays() method is invoked with the new region. This is a recursive 
 mechanism, which seems to display labels in an efficient and balanced way. I 
 am sure that I, or someone else, will make improvements to this mechanism in 
 the future, but it seems to work pretty well.
 
 A companion class exist for CategoryAxisSkipLabels called LabelDefs. This class 
 is constructed by the gatherDrawingGeometries() when the overridden 
 drawCategoryLabels method is called. This method is called whenever a redraw is 
 triggered. The purpose of the LabelDefs class is to hold all the required data 
 (CategoryTick, the Ticks Line marker, the bounds, and anchorPoint) for drawing 
 the labels, and their tick marks. This data is held in a Vector of LabelDef objects. 
 Basically the removeOverlays() method uses this information to mark which labels 
 should be draw, and which should be skipped. This status is marked with a Boolean 
 flag in the LabelDef, which drawLabels checks. Each redraw regenerates the LabelDefs 
 new list. 

 Two flags also compliment the use of CategoryAxisSkipLabels. The first, 
 setTickMarksVisible, is from the super class CategoryAxis. The second flag, 
 displaySkippedTickMarks, is new to this environment. Setting 
 DisplaySkippedTickMarks to true will cause all the tic marks to be drawn, false 
 will draw only the tick marks of drawn labels to be drawn. Of course if I set 
 setTickMarksVisible to false no tick marks are displayed regardless of the 
 DisplaySkippedTickMarks setting.
 ray_lukas@comcast.net || rlukas@pilotsoftware.com

 */
public class CategoryAxisSkipLabels extends CategoryAxis  { 
  
	boolean displaySkippedTickMarks = false;
  
  /*
	 The normal channels for calling this are to truncate the test string. 
	 We wish to avoid this unusual behavior. We do not want this truncating 
	 to happen.. Refer to TextUtilities.createTextBlock() called by 
	 CategoryAxis.createLabel().
  */
  protected TextBlock createLabel(Comparable category,
          RectangleEdge edge, Graphics2D g2) {
	  TextBlock label = createTextBlock(category.toString(), getTickLabelFont(category), 
            getTickLabelPaint(category));
	  return label;
  }
 
  //	we are also changing the way refresh ticks works so that we 
  //	are now calling our 
  public List refreshTicks(Graphics2D g2, AxisState state, Rectangle2D dataArea, RectangleEdge edge) {

	List ticks = new java.util.ArrayList();
	
	// sanity check for data area...
	if (dataArea.getHeight() <= 0.0 || dataArea.getWidth() < 0.0) {
		return ticks;
	}
	
	CategoryPlot plot = (CategoryPlot) getPlot();
	List categories = plot.getCategories();
	double max = 0.0;
	
	if (categories != null) {
		CategoryLabelPosition position = super.getCategoryLabelPositions().getLabelPosition(edge);

		int categoryIndex = 0;
		Iterator iterator = categories.iterator();
		while (iterator.hasNext()) {
			Comparable category = (Comparable) iterator.next();
			
			TextBlock label = createLabel(category, edge, g2);
			if (edge == RectangleEdge.TOP || edge == RectangleEdge.BOTTOM) {
			 max = Math.max(max, 
			         calculateTextBlockHeight(label, position, g2));
			}
			else if (edge == RectangleEdge.LEFT 
			     || edge == RectangleEdge.RIGHT) {
			 max = Math.max(max, 
			         calculateTextBlockWidth(label, position, g2));
			}
			Tick tick = new CategoryTick(category, label, 
			     position.getLabelAnchor(), position.getRotationAnchor(), 
			     position.getAngle());
			ticks.add(tick);
			categoryIndex = categoryIndex + 1;
		}
	}
	state.setMax(max);
	return ticks;
  }

  /*
	  This is where the real skipping work is done. If you want to change the way 
	  that labels are skipped this is the method that you should override. 
	  The basic flow of this method is quite simple. The root function of this 
	  recursive routine is to evaluate the mid point label for a region to see 
	  if it overlaps any other labels that are marked for drawing. If not then 
	  this label is marked for drawing and used in future comparisons. The 
	  current range is split and each subsequent range is pasted into the process 
	  until all ranges have been exhausted.
	
	  Overlapping is discovered by sing a simple intersects between the bounds 
	  of the label we are currently evaluating and its marked for draw neighbors. 
	  Intersects is implemented in the java rectangle class. 
  */
  private LabelDefs removeOverlays(LabelDefs labelDefs, int startTickID, int endTickID){
	  //System.out.println("removeOverlays(labelDefs, " + startTickID + ", " + endTickID + ")");
	  int candidatePosition = ((endTickID - startTickID)/2) + startTickID;
	  //System.out.println("evaluating candidatePosition [" + candidatePosition + "]");
	  Shape candidateShape = labelDefs.getCandidateShape(candidatePosition);
	  if (labelDefs.intersectsAny(candidateShape)) {
		  //System.out.println("intersectAny returned true, break out");
		  return labelDefs;
	  } else {
		  //System.out.println("mark candidatePosition [" + candidatePosition + "] as drawable");
		  labelDefs.markForDraw(candidatePosition);
	  }
	  labelDefs = removeOverlays(labelDefs, startTickID, candidatePosition);	
	  labelDefs = removeOverlays(labelDefs, candidatePosition, endTickID);
	  return labelDefs;
  }
  
  /** 
   * Draws the category labels and returns the updated axis state. 
   * NOTE: This method redefines the corresponding one in <code>CategoryAxis</code>, 
   * and is a copy of that, with added control to skip some labels to be printed. 
   * 
   * @param g2 the graphics device (<code>null</code> not permitted). 
   * @param dataArea the area inside the axes (<code>null</code> not 
   *          permitted). 
   * @param edge the axis location (<code>null</code> not permitted). 
   * @param state the axis state (<code>null</code> not permitted). 
   * @param plotState collects information about the plot (<code>null</code> 
   *          permitted). 
   * 
   * @return The updated axis state (never <code>null</code>). 
   */ 
  protected AxisState drawCategoryLabels(Graphics2D g2, Rectangle2D dataArea, 
                                         RectangleEdge edge, AxisState state, 
                                         PlotRenderingInfo plotState) {
    if (state == null) { 
      throw new IllegalArgumentException("Null 'state' argument."); 
    } 
    
    if (isTickLabelsVisible()) { 
      g2.setFont(getTickLabelFont()); 
      g2.setPaint(getTickLabelPaint()); 
      System.out.println("refresh ticks");
      List ticks = refreshTicks(g2, state, dataArea, edge); 

      try {
    	  state = drawTheChart(g2, dataArea, edge, state, plotState, ticks);
      }
      catch (Exception excep) {
    	  throw new IllegalArgumentException(excep.getMessage()); 
      }
    }
    return state;
  } 
  /*
	  At this point, refreshTicks has been called. We now what to go through 
	  each tick and gather up all the data that we will need to perform the 
	  removeOverlay processing and the eventual drawing. This is precisely 
	  what gatherDrawingGeometries does. 
	  This methods the invokes remove overlays and then performs the label 
	  and tick drawing.
  */
  private AxisState drawTheChart(Graphics2D g2, Rectangle2D dataArea, 
          RectangleEdge edge, AxisState state, 
          PlotRenderingInfo plotState, List ticks) throws Exception {
	  int startTickID = 0;
	  int endTickID = ticks.size();
 
	  LabelDefs labelDefs = gatherDrawingGeometries(g2, dataArea, edge, state, ticks);
	  labelDefs = removeOverlays(labelDefs, startTickID, endTickID);
	  state = drawLabels(g2, edge, state, plotState, dataArea, labelDefs);
	  
	  return state; 
  }
  
 
  // this code derived from Tobi's example
  private Line2D createTickMarkLine(Graphics2D g2, int tickPositon, Rectangle2D dataArea, 
			   RectangleEdge edge, AxisState state){
      float xx = (float) translateValueToJava2D(tickPositon, dataArea, edge); 
      
      Line2D mark = null; 
      double ol = getTickMarkOutsideLength(); 
      double il = getTickMarkInsideLength(); 
      g2.setStroke(getTickMarkStroke()); 
      g2.setPaint(getTickMarkPaint()); 
      if (edge == RectangleEdge.LEFT) { 
          mark = new Line2D.Double(state.getCursor() - ol, xx, state.getCursor() + il, xx); 
      } 
      else if (edge == RectangleEdge.RIGHT) { 
          mark = new Line2D.Double(state.getCursor() + ol, xx, state.getCursor() - il, xx); 
      } 
      else if (edge == RectangleEdge.TOP) { 
          mark = new Line2D.Double(xx, state.getCursor() - ol, xx, state.getCursor() + il); 
      } 
      else if (edge == RectangleEdge.BOTTOM) { 
          mark = new Line2D.Double(xx, state.getCursor() + ol, xx, state.getCursor() - il); 
      } 
      return mark;
  }
  
  private double translateValueToJava2D(int count, 
          Rectangle2D area, 
          RectangleEdge edge) { 
  
     CategoryPlot plot = (CategoryPlot)getPlot(); 
     CategoryAnchor anchor = plot.getDomainGridlinePosition(); 
      RectangleEdge domainAxisEdge = edge;//plot.getDomainAxisEdge(); 
      CategoryDataset data = plot.getDataset(); 
          if (data != null) { 
              CategoryAxis axis = plot.getDomainAxis(); 
              if (axis != null) { 
                  int columnCount = data.getColumnCount(); 
                  return axis.getCategoryJava2DCoordinate( 
                          anchor, count, columnCount, area, domainAxisEdge 
                      ); 
              } 
          } 
  
          return 0.0d; 
  } 

  private Rectangle2D calcAdjustedDataArea(Rectangle2D dataArea, RectangleEdge edge){
      CategoryPlot plot = (CategoryPlot) getPlot(); 
      
      Rectangle2D adjustedDataArea = new Rectangle2D.Double(); 
      if (plot.getRenderer() instanceof Effect3D) { 
          Effect3D e3D = (Effect3D) plot.getRenderer(); 
          double adjustedX = dataArea.getMinX(); 
          double adjustedY = dataArea.getMinY(); 
          double adjustedW = dataArea.getWidth() - e3D.getXOffset(); 
          double adjustedH = dataArea.getHeight() - e3D.getYOffset(); 

          if (edge == RectangleEdge.LEFT || edge == RectangleEdge.BOTTOM) { 
              adjustedY += e3D.getYOffset(); 
          } 
          else if (edge == RectangleEdge.RIGHT || edge == RectangleEdge.TOP) { 
              adjustedX += e3D.getXOffset(); 
          } 
          adjustedDataArea.setRect(adjustedX, adjustedY, adjustedW, adjustedH); 
      } 
      else { 
          adjustedDataArea.setRect(dataArea); 
      } 
      return adjustedDataArea;	  
  }
  private LabelDefs gatherDrawingGeometries(Graphics2D g2, Rectangle2D dataArea, 
          							   RectangleEdge edge, AxisState state, List ticks){
	  int labelCount = 0;
	  LabelDefs labelDefs = new LabelDefs();
	  for (int categoryIndex=0; categoryIndex<ticks.size(); categoryIndex++) {
	      CategoryTick tick = (CategoryTick) ticks.get(categoryIndex); 
	      g2.setPaint(getTickLabelPaint()); 
	
		  CategoryLabelPosition position = getCategoryLabelPositions().getLabelPosition(edge); 
	      Rectangle2D area = calculateNewDrawingRegion(categoryIndex, ticks, dataArea, state, edge);       
	      Point2D anchorPoint = RectangleAnchor.coordinates(area, position.getCategoryAnchor()); 
	
	      TextBlock block = tick.getLabel(); 
	      
	      Shape bounds = block.calculateBounds(g2, (float) anchorPoint.getX(), 
	                                         (float) anchorPoint.getY(), 
	                                         position.getLabelAnchor(), 
	                                         (float) anchorPoint.getX(), 
	                                         (float) anchorPoint.getY(), 
	                                         position.getAngle()); 
	      Line2D mark = createTickMarkLine(g2, categoryIndex, dataArea, edge, state);
	      labelDefs.add("label" + labelCount, tick, mark, bounds, anchorPoint);
	  }
	  return labelDefs;
  }

 
  private AxisState drawLabels(Graphics2D g2, RectangleEdge edge, AxisState state, PlotRenderingInfo plotState, Rectangle2D dataArea, LabelDefs labelDefs)throws Exception {
	  state.setTicks(gatherTicks(labelDefs)); 
	  Object labelDef = null;
	  Iterator labelDefsIter = labelDefs.iterator();
	  
	  while (labelDefsIter.hasNext()) {
		  labelDef = labelDefsIter.next();
		  
		  TextBlock block = ((CategoryTick)labelDefs.getField(labelDef, LabelDefs.CATEGORY_TICK)).getLabel(); 
		  Line2D mark = (Line2D)labelDefs.getField(labelDef, LabelDefs.TICK_MARK_LINE);
		  Point2D anchorPoint = (Point2D)labelDefs.getField(labelDef, LabelDefs.ANCHOR_POINT);
		  CategoryLabelPosition position = getCategoryLabelPositions().getLabelPosition(edge);
		  
		  if ((this.displaySkippedTickMarks) && (super.isTickMarksVisible())){
			  g2.draw(mark);
		  }
		  
		  if (labelDefs.drawThisLabel(labelDef)) {
			  if ((!this.displaySkippedTickMarks) && (super.isTickMarksVisible())){
				  g2.draw(mark);
			  }
		      block.draw(g2, (float) anchorPoint.getX(), (float) anchorPoint.getY(), 
		                 position.getLabelAnchor(), (float) anchorPoint.getX(), 
		                 (float) anchorPoint.getY(), position.getAngle()); 
		  }
	      updatePlotState(plotState, (Shape)labelDefs.getField(labelDef, LabelDefs.BOUNDS));
	  }
	  return updateAxisState(state, edge);
  }


  private List gatherTicks(LabelDefs labelDefs) throws Exception {
	  List ticks = new Vector();
	  Iterator labelDefsIter = labelDefs.iterator();
	  while (labelDefsIter.hasNext()) {
		  ticks.add(labelDefs.getField(labelDefsIter.next(), LabelDefs.CATEGORY_TICK));
	  }
	  return ticks;
  }
  
  private void updatePlotState(PlotRenderingInfo plotState, Shape labelBounds){
		if (plotState != null) { 
			
			ChartRenderingInfo chartRendInfo = plotState.getOwner();
			if (chartRendInfo != null) {
				// seems to be used for processing image maps
				EntityCollection entities = chartRendInfo.getEntityCollection(); 
				if (entities != null) { 
				      //String tooltip = (String) categoryLabelToolTips.get(tick.getCategory()); 
					String tooltip = null; 
					entities.add(new TickLabelEntity(labelBounds, tooltip, null)); 
				} 
			}
		} 
  }
  
  private AxisState updateAxisState(AxisState state, RectangleEdge edge){
      if (edge.equals(RectangleEdge.TOP)) { 
          double h = state.getMax(); 
          state.cursorUp(h); 
        } 
        else if (edge.equals(RectangleEdge.BOTTOM)){ 
          double h = state.getMax(); 
          state.cursorDown(h); 
        } 
        else if (edge == RectangleEdge.LEFT){ 
          double w = state.getMax(); 
          state.cursorLeft(w); 
        } 
        else if (edge == RectangleEdge.RIGHT){ 
          double w = state.getMax(); 
          state.cursorRight(w); 
        } 
      return state;
  }

  
  private Rectangle2D calculateNewDrawingRegion(int categoryIndex, List ticks, Rectangle2D dataArea, 
		  										   AxisState state, RectangleEdge edge) {
      double x0 = 0.0; 
      double x1 = 0.0; 
      double y0 = 0.0; 
      double y1 = 0.0; 
      if (edge == RectangleEdge.TOP) { 
        x0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge); 
        x1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge); 
        y1 = state.getCursor() - getCategoryLabelPositionOffset(); 
        y0 = y1 - state.getMax(); 
      } 
      else if (edge == RectangleEdge.BOTTOM) { 
        x0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge); 
        x1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge); 
        y0 = state.getCursor() + getCategoryLabelPositionOffset(); 
        y1 = y0 + state.getMax(); 
      } 
      else if (edge == RectangleEdge.LEFT) { 
        y0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge); 
        y1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge); 
        x1 = state.getCursor() - getCategoryLabelPositionOffset(); 
        x0 = x1 - state.getMax(); 
      } 
      else if (edge == RectangleEdge.RIGHT) { 
        y0 = getCategoryStart(categoryIndex, ticks.size(), dataArea, edge); 
        y1 = getCategoryEnd(categoryIndex, ticks.size(), dataArea, edge); 
        x0 = state.getCursor() + getCategoryLabelPositionOffset(); 
        x1 = x0 - state.getMax(); 
      } 
      return new Rectangle2D.Double(x0, y0, (x1 - x0), (y1 - y0)); 
  }
  public static TextBlock createTextBlock(final String text, final Font font,
	      final Paint paint) {
	  
		final TextBlock result = new TextBlock();
		result.addLine(text, font, paint);
	
	  return result;
	}
  public void setDisplaySkippedTickMarks(boolean displaySkippedTickMarks){
	  this.displaySkippedTickMarks  = displaySkippedTickMarks;
  }
  
} 