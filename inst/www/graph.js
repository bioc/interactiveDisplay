<style>

.node text {
  pointer-events: none;
  font: 10px sans-serif;
}

.link {
  stroke: #999;
  stroke-opacity: .6;
}

</style>

<script src="http://d3js.org/d3.v2.js"></script>
<script type="text/javascript">var networkOutputBinding = new Shiny.OutputBinding();
  $.extend(networkOutputBinding, {
    find: function(scope) {
      return $(scope).find('.shiny-network-output');
    },
    renderValue: function(el, data) {
      
      //format nodes object
      var nodes = new Array();

      for (var i = 0; i < data.names.length; i++){
        nodes.push({"name":data.names[i],"group":data.groups[i],"title":data.titles[i],"color":data.colors[i]})
      }


      var width = 800;
      var height = 800;
    
	  //var color = d3.scale.category20();
      var lin = data.links
      var force = d3.layout.force()
        .nodes(nodes)
        .links(lin)
        .charge(-120)
        .linkDistance(50)        
        .size([width, height])
        .start();
      
      //remove the old graph
      var svg = d3.select(el).select("svg");      
      svg.remove();
      
      $(el).html("");
      
      //append a new one
      svg = d3.select(el).append("svg");
      
      svg.attr("width", width)
         .attr("height", height);
    
      var link = svg.selectAll("line.link")
          .data(lin)
          .enter().append("line")
          .attr("class", "link")
          .style("stroke-width", function(d) { return Math.sqrt(d.value); });
    
      var node = svg.selectAll("circle.node")
          .data(nodes)
          .enter().append("g")
          .attr("class", "node")
          .call(force.drag);

      node.append("title")
          .attr("dx", 12)
          .attr("dy", ".35em")
          .text(function(d) { return d.title });

      node.append("circle")
          .attr("class", "node")
          .attr("r", 5)
		  //.style("fill", function(d) { return color(d.group); })
      .style("fill", function(d) { return d.color })
		  .on("click", function(d,i) { window.Shiny.onInputChange("probe", d.name);})

	  node.append("text")
		  .attr("dx", 12)
		  .attr("dy", ".35em")
		  .text(function(d) { return d.name });
		        
      force.on("tick", function() {
        link.attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
            .attr("x2", function(d) { return d.target.x; })
            .attr("y2", function(d) { return d.target.y; });
    
        node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
      });
      
    }
  });
  Shiny.outputBindings.register(networkOutputBinding, 'trestletech.networkbinding');
  
</script>

