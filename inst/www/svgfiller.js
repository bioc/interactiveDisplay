<script type="text/javascript">
var networkOutputBinding = new Shiny.OutputBinding();
$.extend(svgOutputBinding, {
    find: function(scope) {
      return $(scope).find('.shiny-network-output');
    },
    renderValue: function(el, data) {
    //$(el).html(data);
      //$(el).html('');
          //use join() to combine the array of strings in data that comprise the svg code
          //data is sent from shiny
          //$(el).html(data.join('')); 
            $(el).html(data);
    //}
    
    //var svg = this.getId(el);
    //var hm = heatmap(el, data);
    
    //var svg = d3.select(el).select.getElementById('gridSVG')
    
    var svg = d3.selectAll("#gridSVG");
    svg.call(d3.behavior.zoom().scaleExtent([1, 8]).on("zoom", zoom));
    
    function zoom() {
      svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
    } 
    }
  });
  Shiny.outputBindings.register(networkOutputBinding, 'balcomes.networkbinding');

</script>
