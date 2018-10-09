// exon class definition
class Exon {
    constructor (graph, column, text = "", length = 1, fullyUTR = false, manX = null, manY = null) {
        // Create a graphic for the exon
        this.graphic = graph.exonLayer.append("ellipse");

        // pointer to graph
        this.graph = graph;

        // Exon position
        this.column = column;
        this.x = column * p.yPadding;
        this.y = p.center;

        // Exon dimensions
        this.radius = p.exonStartSize;
        this.radiusY = p.exonStartSize;
        this.growRadius = this.radius + p.hoverGrow;


        // arrays of exons connected by incoming and outgoing edges
        this.inExons = [];
        this.outExons = [];

        // arrays of exons connected by incoming and outgoing edges
        this.inEdges = [];
        this.outEdges = [];

        // exon properties
        this.textString = text;
        this.length = length;
        this.fullyUTR = fullyUTR;
        this.fill = "white";

        this.manX = manX;
        this.manY = manY;

    }



    // method to render the exon. This is to be called again when the exon properties are updated.
    render() {

        if (this.manX != null) {
            this.x = this.manX;
            this.y = this.manY;
        }

        // set styles if exon is fullyUTR
        if (this.fullyUTR == true) {
            this.graphic.style("stroke-dasharray", ("5, 5"));
            this.graphic.style("fill", "#d9d9d9");
            this.fill = "lightgrey";
        }

        if (this.length < 2) {
            // positioning and scaling settings
            this.graphic.attr("cx", this.x)
                   .attr("cy", this.y)
                   .attr("rx", this.radius)
                   .attr("ry", this.radius)
                   .style("stroke-width", 2)
                   .style("stroke", "black")
                   .style("fill", this.fill);
        } else {
            // if exon spans multiple columns
            // positioning and scaling settings
            this.growRadius = this.radius;
            this.radiusY = (p.padding/2 + 2*p.exonStartSize) / 2;
            this.graphic
                   .attr("cx", this.x + (p.padding/2 + 1*p.exonStartSize) / 2)
                   .attr("cy", this.y)
                   .attr("rx", this.radiusY)
                   .attr("ry", this.radius)
                   .style("stroke-width", 2)
                   .style("stroke", "black")
                   .style("fill", this.fill);
        }

        // get exon data as local variable
        var exonData = this;

        // make graphic grow on hover
        if (this.length == 1) {
            this.graphic.on('mouseover', function(){
                d3.select(this).style("fill", "lightblue");
                d3.select(this).transition().attr("ry", exonData.growRadius).attr("rx", exonData.growRadius).duration(300);
            });
            // make graphic shrink on exit
            this.graphic.on('mouseout', function(){
                d3.select(this).style("fill", exonData.fill);
                d3.select(this).transition().attr("ry", exonData.radius).attr("rx", exonData.radius).duration(300);
            })
        } else {
            this.graphic.on('mouseover', function(){
                d3.select(this).style("fill", "lightblue");
                d3.select(this).transition().attr("ry", exonData.radius + p.hoverGrow).attr("rx", exonData.radiusY + p.hoverGrow).duration(300);
            });
            // make graphic shrink on exit
            this.graphic.on('mouseout', function(){
                d3.select(this).style("fill", exonData.fill);
                d3.select(this).transition().attr("ry", exonData.radius).attr("rx", exonData.radiusY).duration(300);
            })
        }


        // add text
        this.text = this.graph.textLayer.append("text")
           .attr("x", this.graphic.attr("cx"))
           .attr("y", this.graphic.attr("cy"))
           .attr("font-family","sans-serif")
           .attr("pointer-events","none")
           .attr("font-size",12)
           .style("text-anchor", "middle")
           .style("dominant-baseline", "middle")
        //    .text(this.graph.exons.length-1);

        //   if (this.textString.length > 0) {
              this.text.text(this.textString);
            //   this.text.text(this.textString + "," + (this.graph.exons.length-1));
        //   }


        // spread exons at a column out vertically
        for (var column = 1; column < p.numColumns-1; column++) {
            if (this.manX == null) {
                // if there are more than one exons in a column, adjust positions of exons
                var totalHeight = p.yPadding * (this.graph.exonsInColumn[column].length);
                var topExonPos = p.center - (totalHeight * .5);

                // reposition all of the exons at proper positions
                for (var i = 0; i < this.graph.exonsInColumn[column].length; i++) {
                    if (this.graph.exonsInColumn[column][i].manX == null) {
                        this.graph.exonsInColumn[column][i].graphic.attr("cy", topExonPos + (p.yPadding * i));
                        this.graph.exonsInColumn[column][i].text.attr("y", topExonPos + (p.yPadding * i));
                    }
                }
            }
        }
    }

}
