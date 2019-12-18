// on click to get stats
function getStats(){
    var form = $('#statsForm').serializeArray();
    var url = address_base + "/stats/" + form[0]["value"] + "/" + form[1]["value"];
    window.location.replace(url);
}


// var svg = d3.select("#piechart")
// 	.append("svg")
// 	.append("g")

// svg.append("g")
// 	.attr("class", "slices");
// svg.append("g")
// 	.attr("class", "labels");
// svg.append("g")
// 	.attr("class", "lines");

// var width = document.getElementById('piechart').clientWidth,
//     height = width / 3.236,
// 	radius = Math.min(width, height) / 2;

// var pie = d3.layout.pie()
// 	.sort(null)
// 	.value(function(d) {
// 		return d.value;
// 	});

// var arc = d3.svg.arc()
// 	.outerRadius(radius * 0.8)
// 	.innerRadius(radius * 0.4);

// var outerArc = d3.svg.arc()
// 	.innerRadius(radius * 0.9)
// 	.outerRadius(radius * 0.9);

// svg.attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

// var key = function(d){ return d.data.label; };

// var color = d3.scale.ordinal()
// 	.range(["#6b486b", "#d0743c", "#98abc5", "#ff8c00", "#a05d56", "#8a89a6", "#7b6888"]);

// // data should be like this [{label: 'XX', value: '##', color: '#00000'}]    
// var data = [];
// var parsed_info = JSON.parse(info.replace(/\bNaN\b/g, "null"));
// for (var i in parsed_info){  
//         data.push({label: i, value: parsed_info[i][0]})
//     }
// change(data);

// function change(data) {

// 	/* ------- PIE SLICES -------*/
// 	var slice = svg.select(".slices").selectAll("path.slice")
// 		.data(pie(data), key);

// 	slice.enter()
// 		.insert("path")
// 		.style("fill", function(d) { return color(d.data.label); })
// 		.attr("class", "slice");

// 	slice		
// 		.transition().duration(1000)
// 		.attrTween("d", function(d) {
// 			this._current = this._current || d;
// 			var interpolate = d3.interpolate(this._current, d);
// 			this._current = interpolate(0);
// 			return function(t) {
// 				return arc(interpolate(t));
// 			};
// 		})

// 	slice.exit()
// 		.remove();

// 	/* ------- TEXT LABELS -------*/

// 	var text = svg.select(".labels").selectAll("text")
// 		.data(pie(data), key);

// 	text.enter()
// 		.append("text")
// 		.attr("dy", ".35em")
// 		.text(function(d) {
// 			return d.data.label;
// 		});
	
// 	function midAngle(d){
// 		return d.startAngle + (d.endAngle - d.startAngle)/2;
// 	}

// 	text.transition().duration(1000)
// 		.attrTween("transform", function(d) {
// 			this._current = this._current || d;
// 			var interpolate = d3.interpolate(this._current, d);
// 			this._current = interpolate(0);
// 			return function(t) {
// 				var d2 = interpolate(t);
// 				var pos = outerArc.centroid(d2);
// 				pos[0] = radius * (midAngle(d2) < Math.PI ? 1 : -1);
// 				return "translate("+ pos +")";
// 			};
// 		})
// 		.styleTween("text-anchor", function(d){
// 			this._current = this._current || d;
// 			var interpolate = d3.interpolate(this._current, d);
// 			this._current = interpolate(0);
// 			return function(t) {
// 				var d2 = interpolate(t);
// 				return midAngle(d2) < Math.PI ? "start":"end";
// 			};
// 		});

// 	text.exit()
// 		.remove();

// 	/* ------- SLICE TO TEXT POLYLINES -------*/

// 	var polyline = svg.select(".lines").selectAll("polyline")
// 		.data(pie(data), key);
	
// 	polyline.enter()
// 		.append("polyline");

// 	polyline.transition().duration(1000)
// 		.attrTween("points", function(d){
// 			this._current = this._current || d;
// 			var interpolate = d3.interpolate(this._current, d);
// 			this._current = interpolate(0);
// 			return function(t) {
// 				var d2 = interpolate(t);
// 				var pos = outerArc.centroid(d2);
// 				pos[0] = radius * 0.95 * (midAngle(d2) < Math.PI ? 1 : -1);
// 				return [arc.centroid(d2), outerArc.centroid(d2), pos];
// 			};			
// 		});
	
// 	polyline.exit()
// 		.remove();
// };


// data should be like this [{label: 'XX', value: '##', color: '#00000'}]    
var data = [];
var parsed_info = JSON.parse(info.replace(/\bNaN\b/g, "null"));
for (var i in parsed_info){  
        data.push({name: i, value: parsed_info[i][0]})
    }

data.sort(function(a, b){
    return b.value - a.value;
});

var svg = d3.select('svg'),
    canvas = d3.select('#canvas'),
    art = d3.select('#art'),
    labels = d3.select('#labels');

var d3Pie = d3.layout.pie();
d3Pie.value(function(d) {
    return d.value;
});

// store our chart dimensions
var cDim = {
  height: 500,
  width: 500,
  innerRadius: 50,
  outerRadius: 150,
  labelRadius: 180
};

svg.attr({
  height: cDim.height,
  width: '50%'
});

canvas.attr('transform', 'translate(' + (cDim.width/2) + ', ' + (cDim.height/2) + ')');

var pieData = d3Pie(data);

var pieArc = d3.svg.arc()
                .innerRadius(cDim.innerRadius)
                .outerRadius(cDim.outerRadius);


var colors = d3.scale.category20();

var enteringArcs = art.selectAll('.wedge').data(pieData).enter();

enteringArcs.append('path')
            .attr('class', 'wedge')
            .attr('d', pieArc)
            .style('fill', function(d, i){ return colors(i); });



var enteringLabels = labels.selectAll('.label').data(pieData).enter();
var labelGroups = enteringLabels.append('g').attr('class', 'label');

var lines = labelGroups.append('line').attr({
  x1: function(d, i) {
    return pieArc.centroid(d)[0];
  },
  y1: function(d) {
    return pieArc.centroid(d)[1];
  },
  x2: function(d) {
    var centroid = pieArc.centroid(d),
        midAngle = Math.atan2(centroid[1], centroid[0]);
    return Math.cos(midAngle) * cDim.labelRadius;
  },
  y2: function(d) {
    var centroid = pieArc.centroid(d),
        midAngle = Math.atan2(centroid[1], centroid[0]);
    return Math.sin(midAngle) * cDim.labelRadius;
  },
  
  'class': 'label-line',
  'stroke': function(d, i) {
    return colors(i);
  }
});

var textLabels = labelGroups.append('text').attr({
    x: function(d, i) {
      var centroid = pieArc.centroid(d),
          midAngle = Math.atan2(centroid[1], centroid[0]),
          x = Math.cos(midAngle) * cDim.labelRadius,
          sign = x > 0? 1: -1;
      return x + (5*sign);
    },
  
    y: function(d, i) {
      var centroid = pieArc.centroid(d),
          midAngle = Math.atan2(centroid[1], centroid[0]),
          y = Math.sin(midAngle) * cDim.labelRadius;
      return y;
    },
  
    'text-anchor': function(d, i) {
      var centroid = pieArc.centroid(d),
          midAngle = Math.atan2(centroid[1], centroid[0]),
          x = Math.cos(midAngle) * cDim.labelRadius;
      return x > 0? 'start' : 'end';
    },
  
    'class': 'label-text'
}).text(function(d){
  return d.data.name ;
})


// relax the label!
var alpha = 0.5,
    spacing = 15;

function relax() {
  var again = false;
  textLabels.each(function(d, i) {
       var a = this,
           da = d3.select(a),
           y1 = da.attr('y');
        textLabels.each(function(d, j) {
           var b = this;
            if (a === b) {
              return ;
            }
          
            db = d3.select(b);
            if (da.attr('text-anchor') !== db.attr('text-anchor')) {
              return ;
            }
          
            var y2 = db.attr('y');
            deltaY = y1 - y2;
            
            if (Math.abs(deltaY) > spacing) {
              return ;
            }
          
            again = true;
            sign = deltaY > 0? 1: -1;
            var adjust = sign * alpha;
            da.attr('y', +y1 + adjust);
            db.attr('y', +y2 - adjust);
        });
  });
  
  if (again) {
      var labelElements = textLabels[0];
      lines.attr('y2', function(d, i) {
          var labelForLine = d3.select(labelElements[i]);
          return labelForLine.attr('y');
      });
      setTimeout(relax, 20);
  }
}

relax();
