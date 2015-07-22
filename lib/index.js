
var parser = require("biojs-io-xmfa");
var fastaParser = require("biojs-io-fasta");
var brewer = require("colorbrewer");
var d3 = require("d3");

var biojsvizxmfa;

module.exports = biojsvizxmfa = function(opts){
  this.el = opts.el;
  this.opts = {
    "x_offset": 10,
    "y_offset": 10,
    "scale": 4,
    "genome_sep": 30,
    "color": brewer['Set1'][8]
  };
  var self = this;

  parser.read("/outseq.xmfa", function(err, model) {
    var text = document.createElement("div");

    // Convert lcbs to easy-to-read table
    var table = document.createElement("div");
    table.innerHTML = "<h1>LCB Table</h1>" + biojsvizxmfa.lcbsAsTable(self, model);
    text.appendChild(table);

    var svgContainer = d3.select(self.el).append("svg")
                                         .attr("width", 1000)
                                         .attr("height", 200);

    fastaParser.read("/outseq.fa", function(err, fa_model) {
      var lcb_fa_table = {};
      for(var fa_idx in fa_model){
        var seq = fa_model[fa_idx];
        var seq_seq = seq.seq;

        lcb_fa_table[1 + parseInt(fa_idx)] = seq.name;

        var line_attr = {
          "x1": self.opts.x_offset + 0,
          "y1": self.opts.y_offset + self.opts.genome_sep * fa_idx,
          "x2": self.opts.x_offset + (seq_seq.length / self.opts.scale),
          "y2": self.opts.y_offset + self.opts.genome_sep * fa_idx,
          "stroke-width": 2,
          "stroke": "black"
        };
        var line_label_attr = {
          "dx": 0,
          "dy": "-.35em",
          "font-family": "sans-serif",
          "font-size": "10px",
          "fill": "black"
        };
        var line = svgContainer.append("line")
                               .attr(line_attr)
                               .append("text")
                               .attr("dx", 12)
                               .attr("dy", ".35em")
                               .text("Blah");
                               //.text(seq.name)
                               //.attr(line_label_attr);

      }

      for(var lcb_idx in model){if(model[lcb_idx].length > 1){ 
        var lcb_group = svgContainer.append("g");

        // Find a nice colour
        var current_colour = self.opts.color[lcb_idx];

        // Store polygon points as left/right sides
        var current_poly_left = [];
        var current_poly_right = [];
        var regions = [];

        // Loop across Regions in LCB
        // Plot them + record left/right sides
        for(var lcb_subidx in model[lcb_idx]){
          var hsp = model[lcb_idx][lcb_subidx];
          var region_highlight_attr = {
            "x": self.opts.x_offset + (hsp['start'] / self.opts.scale),
            "width": ((hsp['end'] - hsp['start']) / self.opts.scale),
            "y": self.opts.y_offset + self.opts.genome_sep * (hsp['lcb_idx'] - 1) - 3,
            "height": 6,
            "fill": current_colour
          }
          regions.push(region_highlight_attr);

          var left_x, right_x;
          if(hsp['strand'] == '+'){
            left_x = hsp['start'];
            right_x = hsp['end'];
          }else{
            left_x = hsp['end'];
            right_x = hsp['start'];
          }

          current_poly_left.push({
            'x': self.opts.x_offset + (left_x / self.opts.scale),
            'y': self.opts.y_offset + self.opts.genome_sep * (hsp['lcb_idx'] - 1)
          });
          current_poly_right.push({
            'x': self.opts.x_offset + (right_x / self.opts.scale),
            'y': self.opts.y_offset + self.opts.genome_sep * (hsp['lcb_idx'] - 1)
          });
        }

        current_poly_left.reverse();
        for(var tmp in current_poly_left){
          current_poly_right.push(current_poly_left[tmp]);
        }

        var poly = lcb_group.append("polygon")
                            .data([current_poly_right])
                            .attr("points",function(d) { 
                               return d.map(function(d) {
                                 return [d.x, d.y].join(",");
                               }).join(" ");
                            })
                            .attr("opacity", 0.6)
                            .attr("fill", current_colour);


        for(var idx in regions){
          lcb_group.append("rect").attr(regions[idx]);
        }


      }}

      svgContainer.selectAll("g")
                  .on("mouseover", function(){
                     d3.select(this).select("polygon")
                       .style("opacity", 0.6)
                       .transition()
                       .duration(500)
                       .style("opacity", 1.0);
 
                     d3.select(this).selectAll("rect")
                       .style("stroke-width", 0)
                       .style("stroke", "black")
                       .transition()
                       .duration(500)
                       .style("stroke-width", 1)
                       .style("stroke", "black");

                     var element = d3.select(this)[0][0];
                     element.parentNode.appendChild(element);

                   })
                   .on("mouseout", function() {
                     d3.select(this).select("polygon")
                       .style("opacity", 1.0)
                       .transition()
                       .duration(500)
                       .style("opacity", 0.6);

                     d3.select(this).selectAll("rect")
                       .style("stroke-width", 1)
                       .style("stroke", "black")
                       .transition()
                       .duration(500)
                       .style("stroke-width", 0)
                       .style("stroke", "black");
                   });



    });
    self.el.appendChild(text);
  });
};


biojsvizxmfa.lcbsAsTable = function (self, lcbs) {
  var table_text = "";
  
  for(var lcb_idx in lcbs){if(lcbs[lcb_idx].length > 1){
    for(var lcb_subidx in lcbs[lcb_idx]){
      var hsp = lcbs[lcb_idx][lcb_subidx];
      var row = [lcb_idx, hsp['lcb_idx'], hsp['start'], hsp['end'], hsp['strand']];
      table_text += "<tr style='background:" + self.opts.color[lcb_idx] + "'><td>" + row.join("</td><td>") + "</td></tr>";
    }
  }}
  var thead = "<thead><tr><th>" + ["LCB", "Sequence", "Start", "End", "Strand"].join("</th><th>") + " </th></tr></thead>";
  return "<table>" + thead + "<tbody>" + table_text + "</tbody></table>";
};

