var region = '';
var regionlist = [];
var data = {};


function entry(region, regionlist, passedregions, regioncoords, sequences, ctx) {


    window.plot = null;

    window.normtotal = false;
    window.normmax = false;

    window.region = region;
    window.regionlist = regionlist;
    window.passedregions = passedregions;
    window.regioncoords = regioncoords;
    window.sequences = sequences;
    window.ctx = ctx;

    window.y_profile = 'COV';
    window.y2_profile = null;

    $.jqplot.config.enablePlugins = true;

    readProfilesData();

    makePlot();

    window.plot.target.bind('jqplotZoom', function(ev, gridpos, datapos, plot, cursor){
        handleZooming();
    });

    loadSideList();

    var rowToSelect = $("td").filter(function() {
        return $(this).text() == region;
    }).closest("tr");

    swapRowSelection(rowToSelect);

    var container = $('#scroller');
    container.scrollTop(0);
    container.scrollTop( rowToSelect.offset().top - container.offset().top - 220);


};

function makePlot(){

    data1 = prepareDataForPlotting(window.data[window.y_profile]);

    $("#title").html(window.region+' ('+regionAsString(window.regioncoords[window.region])+')');

    window.plot = $.jqplot('myplot', [data1],
        {

                cursor: {
					show: true,
					zoom:true,
					looseZoom: true,
					showTooltip:false,
					showTooltipUnitPosition:true,
					constrainZoomTo:'x',
					constrainOutsideZoom:true,
					showVerticalLine:true,
					dblClickReset:false,
					tooltipLocation:'sw'
				},

                highlighter: {
					show: true,
					sizeAdjust: 12,
					tooltipAxes: "xy",
					tooltipLocation: 'n',
					formatString: "Position: %d; Value: %g",
					bringSeriesToFront:true,
                    fadeTooltip: true,
                    tooltipFadeSpeed: "fast",
				},

                seriesDefaults: {
                    showMarker: true,
                    pointLabels: { show: false }
                },

                axes: {
					xaxis: {
						min: window.regioncoords[window.region][1],
				 		max: window.regioncoords[window.region][2],
                        numberTicks: 6,
                        tickOptions: { showLabel:true, formatString: '%d', fontSize: 14 }
					},
                    yaxis: {
						min: 0,
                        max: getMaxY(data1)*1.1,
                        numberTicks: 6,
                        tickOptions: { showLabel:true, formatString: "%.2f", fontSize: 14 }
					},
					y2axis: {
						min: 0
					}
				}
        }
    );

    //window.plot.redraw();

    window.plot.replot();

    var x=$("#myplot").offset().left;
    var y=$("#myplot").offset().top;
    $('#canvas').offset({top:y+500,left:x});




};

function regionAsString(coords){
    return coords[0]+':'+coords[1].toString()+'-'+coords[2].toString();
};


function handleZooming(){
    var plotData =  window.plot.series[0].data;
    var maxy = -1;
    var miny = 9999999;
    var totalmaxy = -1;
    for (var i=0; i<plotData.length; i++){
        if (plotData[i][1]>totalmaxy)
                totalmaxy=plotData[i][1];
        if (plotData[i][0] >= window.plot.axes.xaxis.min && plotData[i][0] <= window.plot.axes.xaxis.max){
            if (plotData[i][1]>maxy)
                maxy=plotData[i][1];
            if (plotData[i][1]<miny)
                miny=plotData[i][1];
        }
    }

    if (maxy==miny){
        var w=maxy;
        miny=w*0.8;
        maxy=w*1.2;
    }

    var maxy2=-1;
    var miny2=9999999;
    var totalmaxy2 = -1;
    if (window.plot.series.length==2) {
        var plotData2 =  window.plot.series[1].data;
        for (var i=0; i<plotData.length; i++){
            if (plotData2[i][1]>totalmaxy2)
                totalmaxy2=plotData2[i][1];
            if (plotData[i][0] >= window.plot.axes.xaxis.min && plotData[i][0] <= window.plot.axes.xaxis.max){
                if (plotData2[i][1]>maxy2)
                    maxy2=plotData2[i][1];
                if (plotData2[i][1]<miny2)
                    miny2=plotData2[i][1];
            }
        }
    }

    if (maxy2==miny2){
        var w=maxy2;
        miny2=w*0.8;
        maxy2=w*1.1;
    }


    scaleXAxis();

    scaleYAxes(miny,maxy,miny2,maxy2,totalmaxy,totalmaxy2);

    makeReferenceBar(window.sequences[window.region]);

    $("#zoombuttons").show();
};


function get_x_array(){
    var plotData =  window.plot.series[0].data;
    var ret = [];
    for (var i=0; i<plotData.length; i++){
        if (plotData[i][0] >= window.plot.axes.xaxis.min && plotData[i][0] <= window.plot.axes.xaxis.max)
            ret.push(plotData[i][0]);
    }
    return ret
}


function scaleXAxis(){

    x_array = get_x_array();

    var first = x_array[0];
    var last = x_array[x_array.length-1];

    var interval=Math.floor((last-first)/5);
    if ((last-first)>=10){
        window.plot.axes.xaxis.tickInterval=interval;
        window.plot.axes.xaxis.numberTicks=Math.floor((last-first)/interval)+1;
    }

    if ((last-first)<10){
        window.plot.axes.xaxis.min=first;
        window.plot.axes.xaxis.max=last;
        window.plot.axes.xaxis.tickInterval=1;
        window.plot.axes.xaxis.numberTicks=last-first+1;
    }
};



function makeReferenceBar(sequence) {

    x_array = get_x_array();

    var width = document.getElementById('myplot').offsetWidth + 30;
    window.ctx.canvas.width = width;


    window.ctx.clearRect(0,0,width,20);
    //window.ctx.fillStyle='#E0E0E0';
    //window.ctx.fillRect(0,0,width,20);

    var regioncoords = window.regioncoords[window.region];

    for (var i=0; i<x_array.length; i++) {

        var base=sequence[x_array[i]-regioncoords[1]];


  			   window.ctx.fillStyle='#C0C0C0';
  			   if (base=="A"){
				   window.ctx.fillStyle='green';
  			   }
			   if (base=="C"){
				   window.ctx.fillStyle='blue';
			   }
			   if (base=="G"){
				   window.ctx.fillStyle='black';
			   }
			   if (base=="T"){
				   window.ctx.fillStyle='red';
			   }

				var dx=window.plot.axes.xaxis.u2p(x_array[i]);

				if (x_array.length==1) {
					w=1000;
				}

				if (i>0) {
					var dxprev=window.plot.axes.xaxis.u2p(x_array[i-1]);
					w=(dx-dxprev)/2;
				}

				if (i<x_array.length-1) {
					var dxnext=window.plot.axes.xaxis.u2p(x_array[i+1]);
					w=(dxnext-dx)/2;
				}

				window.ctx.fillRect(dx-w,0,2*w,20);

				if (x_array.length<=40){
					window.ctx.lineWidth = 0.5;
					window.ctx.strokeStyle = 'white';
					window.ctx.strokeRect(dx-w,0,2*w,20);

					window.ctx.fillStyle='white';
					window.ctx.font = "16px Arial";
					if (base=="A"){
						window.ctx.fillText("A",dx-6,17);
					}
					if (base=="C"){
						window.ctx.fillText("C",dx-6,17);
					}
					if (base=="G"){
						window.ctx.fillText("G",dx-6,17);
					}
					if (base=="T"){
						window.ctx.fillText("T",dx-6,17);
					}
				}
    }

};


function scaleYAxes(miny1,maxy1,miny2,maxy2,totalmaxy,totalmaxy2){

    if (window.normtotal){
        window.plot.axes.yaxis.max=totalmaxy;
        window.plot.axes.yaxis.min=0;
        window.plot.axes.y2axis.max=totalmaxy2;
        window.plot.axes.y2axis.min=0;
    }
    else {
        window.plot.axes.yaxis.max=maxy1;
        window.plot.axes.yaxis.min=miny1;
        window.plot.axes.y2axis.max=maxy2;
        window.plot.axes.y2axis.min=miny2;
    }

    if (window.normmax){
        window.plot.axes.yaxis.max=Math.max(window.plot.axes.yaxis.max,window.plot.axes.y2axis.max)
        window.plot.axes.y2axis.max=window.plot.axes.yaxis.max
        window.plot.axes.yaxis.min=Math.min(window.plot.axes.yaxis.min,window.plot.axes.y2axis.min)
        window.plot.axes.y2axis.min=window.plot.axes.yaxis.min
    }


    var delta=(window.plot.axes.yaxis.max-window.plot.axes.yaxis.min)*0.1;
    window.plot.axes.yaxis.min=Math.max(window.plot.axes.yaxis.min-delta,0);
    window.plot.axes.yaxis.max=window.plot.axes.yaxis.max+delta;
    window.plot.axes.yaxis.tickInterval=(window.plot.axes.yaxis.max-window.plot.axes.yaxis.min)/5;
    window.plot.axes.yaxis.numberTicks=6;

    var delta=(window.plot.axes.y2axis.max-window.plot.axes.y2axis.min)*0.1;
    window.plot.axes.y2axis.min=Math.max(window.plot.axes.y2axis.min-delta,0);
    window.plot.axes.y2axis.max=window.plot.axes.y2axis.max+delta;
    window.plot.axes.y2axis.tickInterval=(window.plot.axes.y2axis.max-window.plot.axes.y2axis.min)/5;
    window.plot.axes.y2axis.numberTicks=6;

    window.plot.replot();

};


function loadSideList(){

    var rows = '';

    $.each(window.regionlist, function(index, item) {
        var row = '';

        if (window.passedregions.indexOf(item) > -1)
            row += '<tr onclick="swapRowSelection(this);">';
        else
            row += '<tr onclick="swapRowSelection(this);" style="color: #e22c11; font-weight: bold;">';

        row += '<td class="column1">' + item + '</td>';
        rows += row + '</tr>';

    });
    $("#contentrows").html(rows);

};


function swapRowSelection(row) {
    $("#contentrows tr").removeClass("selected");
    $(row).addClass("selected");
    window.region = whichSelected();
    readProfilesData();
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
};


function prepareDataForPlotting(profile){

    var start = window.regioncoords[window.region][1];
    var end = window.regioncoords[window.region][2];

    var ret = [];
    for (var i = start; i < end; i++) {
        ret.push([i,profile[i-start]]);
    };

    return ret;
};


function whichSelected(){
    var sel = '';
    $('.selected').each(function() {
        sel = $(this).find('td:first').text();
    });
    return sel
};


function getMaxY(data){
    var maxy=-1;
    for (var i = 0; i < data.length; i++) {
        var y=data[i][1];
        if (y>maxy){
            maxy=y;
        }
    }
    return maxy;
};


function readProfilesData(){
    $.ajax({
            type : 'POST',
            url : "/profiles",
            contentType: 'application/json;charset=UTF-8',
            data : JSON.stringify({ region: window.region }),
            success: function(data) {
               window.data = data;
            },
            async: false
    });
};

function handleZoomOut(){
    window.plot.resetZoom();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
};


function goLeft(){

    var plotData =  window.plot.series[0].data;
    if (window.plot.axes.xaxis.min - 1 < plotData[0][0])
        return;

    window.plot.axes.xaxis.min=window.plot.axes.xaxis.min-1;
    window.plot.axes.xaxis.max=window.plot.axes.xaxis.max-1;
    window.plot.replot();
    makeReferenceBar(window.sequences[window.region]);
};

function goRight(){

    var plotData =  window.plot.series[0].data;
    if (window.plot.axes.xaxis.max + 1 > plotData[plotData.length-1][0])
        return;

    window.plot.axes.xaxis.min=window.plot.axes.xaxis.min+1;
    window.plot.axes.xaxis.max=window.plot.axes.xaxis.max+1;
    window.plot.replot();
    makeReferenceBar(window.sequences[window.region]);
};

function changeMetrics(){
    window.y_profile = $("#metrics :selected").text();
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
};

