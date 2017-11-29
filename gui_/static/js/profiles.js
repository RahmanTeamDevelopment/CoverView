var region = '';
var regionlist = [];
var data = {};


function entry(region, regionlist, passedregions, regioncoords, sequences, ctx, passdef) {


    window.plot = null;

    window.norm = true;

    window.reads = 'all';

    if (region=='')
        region = regionlist[0];

    window.region = region;
    window.regionlist = regionlist;


    window.passedregions = passedregions;
    window.regioncoords = regioncoords;
    window.sequences = sequences;
    window.ctx = ctx;
    window.passdef = passdef;
    window.showoverlay = true;

    window.yaxis_min_saved = null;
    window.yaxis_max_saved = null;
    window.yaxis_tickinterval_saved = null;
    window.yaxis_numberticks_saved = null;

    window.y2axis_min_saved = null;
    window.y2axis_max_saved = null;
    window.y2axis_tickinterval_saved = null;
    window.y2axis_numberticks_saved = null;

    window.y_profile = 'COV';
    window.y2_profile = '---';

    window.zoomedin = false;

    $.jqplot.config.enablePlugins = true;

    readProfilesData();


    if (!('COV+' in window.data)){
        document.getElementById("allreads").disabled=true;
        $("#forwardlabel").css("color", "#a4abb5");
        document.getElementById("forward").disabled=true;
        $("#reverselabel").css("color", "#a4abb5");
        document.getElementById("reverse").disabled=true;
    }



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
    container.scrollTop( rowToSelect.offset().top - container.offset().top - 200);


};


function singleGraphPlot(data1, cutoff_y){

    var showcutoff = true;
    var cutofftxt = '';
    if (cutoff_y == '.')
        showcutoff = false;
    else
        cutofftxt = window.y_profile+'='+cutoff_y.toString();


    var overlay = {};
    if (window.showoverlay && window.reads=='all') {
        overlay = {
            show: true,
            objects: [
                {
                    dashedHorizontalLine: {
                        show: showcutoff,
                        name: 'cut1line',
                        y: cutoff_y,
                        lineWidth: 3,
                        color: '#00749F',
                        shadow: false,
                        xOffset: 0,
                        showTooltip: showcutoff,
                        tooltipFormatString: cutofftxt
                    }
                }
            ]
        };
    }


    return $.jqplot('myplot', [data1],
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
					formatString: "%d; %g",
					bringSeriesToFront:true,
                    fadeTooltip: true,
                    tooltipFadeSpeed: "fast",
				},

                seriesColors: ['#00749F'],

                canvasOverlay: overlay,

                seriesDefaults: {
                    showMarker: false,
                    pointLabels: { show: false },
                    lineWidth:3
                },

                axes: {
					xaxis: {
						min: window.regioncoords[window.region][1],
				 		max: window.regioncoords[window.region][2]-1,
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
}


function doubleGraphPlot(data1, data2, cutoff_y, cutoff_y2){

    var objs = [];
    if (cutoff_y != '.'){
        objs.push({
                    dashedHorizontalLine: {
                        show: true,
                        name: 'cut1line',
                        y: cutoff_y,
                        lineWidth: 3,
                        color: '#00749F',
                        shadow: false,
                        xOffset: 0,
                        showTooltip: true,
                        tooltipFormatString: window.y_profile+'='+cutoff_y.toString()
                    }
                });
    }

    if (cutoff_y2 != '.'){
        objs.push({
                    dashedHorizontalLine: {
                        show: true,
                        name: 'cut2line',
                        y: cutoff_y2,
                        yaxis: 'y2axis',
                        lineWidth: 3,
                        color: 'darkred',
                        shadow: false,
                        xOffset: 0,
                        showTooltip: true,
                        tooltipFormatString: window.y2_profile+'='+cutoff_y2.toString()
                    }
                });
    }

    var overlay = {};
    if (window.showoverlay && window.reads=='all') {
        overlay = {
            show: true,
            objects: objs
        };
    }


    return $.jqplot('myplot', [data1, data2],
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

                seriesColors: ['#00749F', 'darkred'],

                highlighter: {
					show: true,
					sizeAdjust: 12,
					tooltipAxes: "xy",
					tooltipLocation: 'n',
					bringSeriesToFront:true,
                    fadeTooltip: true,
                    tooltipFadeSpeed: "fast",
                    formatString: "%d; %g",
				},

                canvasOverlay: overlay,

                seriesDefaults: {
                    showMarker: false,
                    pointLabels: { show: false },
                    lineWidth:3
                },

                axes: {
					xaxis: {
						min: window.regioncoords[window.region][1],
				 		max: window.regioncoords[window.region][2]-1,
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
						min: 0,
                        max: getMaxY(data2)*1.1,
                        numberTicks: 6,
                        tickOptions: { showLabel:true, formatString: "%.2f", fontSize: 14 }
					}
				},

                series:
                    [
						{yaxis:'yaxis'},
						{yaxis:'y2axis'}
                    ]
        }
    );
}

function makePlot(){


    var cutoff_y = '.';

    if (window.y_profile=='COV'){
        if ('MEDCOV_MIN' in window.passdef)
            cutoff_y = Number(window.passdef['MEDCOV_MIN']);
        if ('MINCOV_MIN' in window.passdef)
            cutoff_y = Number(window.passdef['MINCOV_MIN']);
    }

    if (window.y_profile=='QCOV'){
        if ('MEDQCOV_MIN' in window.passdef)
            cutoff_y = Number(window.passdef['MEDQCOV_MIN']);
        if ('MINQCOV_MIN' in window.passdef)
            cutoff_y = Number(window.passdef['MINQCOV_MIN']);
    }

    if (window.y_profile=='FLBQ'){
        if ('MAXFLBQ_MAX' in window.passdef)
            cutoff_y = Number(window.passdef['MAXFLBQ_MAX']);
    }

    if (window.y_profile=='FLMQ'){
        if ('MAXFLMQ_MAX' in window.passdef)
            cutoff_y = Number(window.passdef['MAXFLMQ_MAX']);
    }

    var cutoff_y2 = '.';

    if (window.y2_profile=='COV'){
        if ('MEDCOV_MIN' in window.passdef)
            cutoff_y2 = Number(window.passdef['MEDCOV_MIN']);
        if ('MINCOV_MIN' in window.passdef)
            cutoff_y2 = Number(window.passdef['MINCOV_MIN']);
    }

    if (window.y2_profile=='QCOV'){
        if ('MEDQCOV_MIN' in window.passdef)
            cutoff_y2 = Number(window.passdef['MEDQCOV_MIN']);
        if ('MINQCOV_MIN' in window.passdef)
            cutoff_y2 = Number(window.passdef['MINQCOV_MIN']);
    }

    if (window.y2_profile=='FLBQ'){
        if ('MAXFLBQ_MAX' in window.passdef)
            cutoff_y2 = Number(window.passdef['MAXFLBQ_MAX']);
    }

    if (window.y2_profile=='FLMQ'){
        if ('MAXFLMQ_MAX' in window.passdef)
            cutoff_y2 = Number(window.passdef['MAXFLMQ_MAX']);
    }

    var prof = window.y_profile;
    if (window.reads == 'forward'){
        prof = prof + '+';
    }
    if (window.reads == 'reverse'){
        prof = prof + '-';
    }

    var prof2 = window.y2_profile;
    if (window.reads == 'forward'){
        prof2 = prof2 + '+';
    }
    if (window.reads == 'reverse'){
        prof2 = prof2 + '-';
    }


    data1 = prepareDataForPlotting(window.data[prof]);

    if (window.y2_profile == '---'){
        window.plot = singleGraphPlot(data1, cutoff_y);
        if (window.plot.axes.yaxis.max ==0 && window.plot.axes.yaxis.min ==0){
            window.plot.axes.yaxis.min=0;
            window.plot.axes.yaxis.max=1;
        }
    }
    else {
        data2 = prepareDataForPlotting(window.data[prof2]);
        window.plot = doubleGraphPlot(data1, data2, cutoff_y, cutoff_y2);
        if (((window.y_profile == 'COV') && (window.y2_profile == 'QCOV')) || ((window.y_profile == 'QCOV') && (window.y2_profile == 'COV'))) {
            var m = Math.max(window.plot.axes.yaxis.max, window.plot.axes.y2axis.max);
            window.plot.axes.yaxis.max = m;
            window.plot.axes.y2axis.max = m;
            window.plot.axes.yaxis.tickInterval = (window.plot.axes.yaxis.max - window.plot.axes.yaxis.min) / 5;
            window.plot.axes.y2axis.tickInterval = (window.plot.axes.y2axis.max - window.plot.axes.y2axis.min) / 5;
        }

        if (window.plot.axes.yaxis.max ==0 && window.plot.axes.yaxis.min ==0){
            window.plot.axes.yaxis.min=0;
            window.plot.axes.yaxis.max=1;
        }
        if (window.plot.axes.y2axis.max ==0 && window.plot.axes.y2axis.min ==0){
            window.plot.axes.y2axis.min=0;
            window.plot.axes.y2axis.max=1;
        }

    }



    $("#title").html(window.region+' ('+regionAsString(window.regioncoords[window.region])+')');

    window.plot.replot();

    var x=$("#myplot").offset().left;
    var y=$("#myplot").offset().top;
    $('#canvas').offset({top:y+500,left:x});

    save_normal_scale_values();


};

function regionAsString(coords){
    return coords[0]+':'+coords[1].toString()+'-'+coords[2].toString();
};


function handleZooming(){

    window.zoomedin = true;

    if (window.norm) {

        normalizeY();
        if (window.y2_profile != '---') {
            normalizeY2();
        }

    }


    scaleXAxis();

    window.plot.replot();

    makeReferenceBar(window.sequences[window.region]);

    x_array = get_x_array();
    $("#zoomedcoords").text(x_array[0]+"-"+x_array[x_array.length-1]);

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

    window.plot.axes.xaxis.min = first;
    window.plot.axes.xaxis.max = last;

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

    if (window.plot.axes.yaxis.max ==0 && window.plot.axes.yaxis.min ==0){
            window.plot.axes.yaxis.min=0;
            window.plot.axes.yaxis.max=1;
    }
    if (window.y2_profile != '---'){
        if (window.plot.axes.y2axis.max ==0 && window.plot.axes.y2axis.min ==0){
            window.plot.axes.y2axis.min=0;
            window.plot.axes.y2axis.max=1;
        }
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
    window.zoomedin = false;
};


function prepareDataForPlotting(profile){

    var start = window.regioncoords[window.region][1];
    var end = window.regioncoords[window.region][2];

    var ret = [];
    for (var i = start; i < end; i++) {

        var v = profile[i-start];

        if (v==null)
            v = 0;

        ret.push([i,v]);
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

     if (window.plot.axes.yaxis.max ==0 && window.plot.axes.yaxis.min ==0){
            window.plot.axes.yaxis.min=0;
            window.plot.axes.yaxis.max=1;
    }
    if (window.y2_profile != '---'){
        if (window.plot.axes.y2axis.max ==0 && window.plot.axes.y2axis.min ==0){
            window.plot.axes.y2axis.min=0;
            window.plot.axes.y2axis.max=1;
        }
    }
    window.plot.replot();

    if (((window.y_profile == 'COV') && (window.y2_profile == 'QCOV')) || ((window.y_profile == 'QCOV') && (window.y2_profile == 'COV'))) {
        var m = Math.max(window.plot.axes.yaxis.max, window.plot.axes.y2axis.max);
        window.plot.axes.yaxis.max = m;
        window.plot.axes.y2axis.max = m;
        window.plot.axes.yaxis.tickInterval = (window.plot.axes.yaxis.max - window.plot.axes.yaxis.min) / 5;
        window.plot.axes.y2axis.tickInterval = (window.plot.axes.y2axis.max - window.plot.axes.y2axis.min) / 5;
        window.plot.replot();
    }

    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();

    window.zoomedin = false;
};


function goLeft(){

    var plotData =  window.plot.series[0].data;
    if (window.plot.axes.xaxis.min - 1 < plotData[0][0])
        return;

    window.plot.axes.xaxis.min=window.plot.axes.xaxis.min-1;
    window.plot.axes.xaxis.max=window.plot.axes.xaxis.max-1;

    if (window.norm) {

        normalizeY();
        if (window.y2_profile != '---') {
            normalizeY2();
        }

    }

    if (window.plot.axes.yaxis.max ==0 && window.plot.axes.yaxis.min ==0){
            window.plot.axes.yaxis.min=0;
            window.plot.axes.yaxis.max=1;
    }
    if (window.y2_profile != '---'){
        if (window.plot.axes.y2axis.max ==0 && window.plot.axes.y2axis.min ==0){
            window.plot.axes.y2axis.min=0;
            window.plot.axes.y2axis.max=1;
        }
    }

    window.plot.replot();
    makeReferenceBar(window.sequences[window.region]);

    x_array = get_x_array();
    $("#zoomedcoords").text(x_array[0]+"-"+x_array[x_array.length-1]);
};

function goRight(){

    var plotData =  window.plot.series[0].data;
    if (window.plot.axes.xaxis.max + 1 > plotData[plotData.length-1][0])
        return;

    window.plot.axes.xaxis.min=window.plot.axes.xaxis.min+1;
    window.plot.axes.xaxis.max=window.plot.axes.xaxis.max+1;

    if (window.norm) {
        normalizeY();
        if (window.y2_profile != '---') {
            normalizeY2();
        }
    }

    if (window.plot.axes.yaxis.max ==0 && window.plot.axes.yaxis.min ==0){
            window.plot.axes.yaxis.min=0;
            window.plot.axes.yaxis.max=1;
    }
    if (window.y2_profile != '---'){
        if (window.plot.axes.y2axis.max ==0 && window.plot.axes.y2axis.min ==0){
            window.plot.axes.y2axis.min=0;
            window.plot.axes.y2axis.max=1;
        }
    }


    window.plot.replot();
    makeReferenceBar(window.sequences[window.region]);

    x_array = get_x_array();
    $("#zoomedcoords").text(x_array[0]+"-"+x_array[x_array.length-1]);
};

function changeMetrics(){
    window.y_profile = $("#metrics :selected").text();
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
    window.zoomedin = false;
};

function changeMetrics2(){
    window.y2_profile = $("#metrics2 :selected").text();
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
    window.zoomedin = false;
};

function normalizeY(){

    if (!window.zoomedin)
        return '';

    var plotData =  window.plot.series[0].data;
    var y=[];
    for (i=0;i<plotData.length;i++) {
        if (plotData[i][0] >= window.plot.axes.xaxis.min && plotData[i][0] <= window.plot.axes.xaxis.max)
            y.push(plotData[i][1]);
    }
    var maxtotal_y = Math.max.apply(null, y);
    var mintotal_y = Math.min.apply(null, y);

    if (maxtotal_y==mintotal_y){
        var w=maxtotal_y;
        mintotal_y=w*0.8;
        maxtotal_y=w*1.2;
    }

    window.plot.axes.yaxis.max=maxtotal_y;
    window.plot.axes.yaxis.min=mintotal_y;

    var delta=(window.plot.axes.yaxis.max-window.plot.axes.yaxis.min)*0.1;
    window.plot.axes.yaxis.min=Math.max(window.plot.axes.yaxis.min-delta,0);
    window.plot.axes.yaxis.max=window.plot.axes.yaxis.max+delta;
    window.plot.axes.yaxis.tickInterval=(window.plot.axes.yaxis.max-window.plot.axes.yaxis.min)/5;
    window.plot.axes.yaxis.numberTicks=6;

     if (window.plot.axes.yaxis.max ==0 && window.plot.axes.yaxis.min ==0){
            window.plot.axes.yaxis.min=0;
            window.plot.axes.yaxis.max=1;
    }


};


function save_normal_scale_values(){

    window.yaxis_min_saved = window.plot.axes.yaxis.min;
    window.yaxis_max_saved = window.plot.axes.yaxis.max;
    window.yaxis_tickinterval_saved = window.plot.axes.yaxis.tickInterval;
    window.yaxis_numberticks_saved = window.plot.axes.yaxis.numberTicks;

    if (window.y2_profile != '---') {
        window.y2axis_min_saved = window.plot.axes.y2axis.min;
        window.y2axis_max_saved = window.plot.axes.y2axis.max;
        window.y2axis_tickinterval_saved = window.plot.axes.y2axis.tickInterval;
        window.y2axis_numberticks_saved = window.plot.axes.y2axis.numberTicks;
    }

};


function normal_scaleY(){
    window.plot.axes.yaxis.min = window.yaxis_min_saved;
    window.plot.axes.yaxis.max = window.yaxis_max_saved;
    window.plot.axes.yaxis.tickInterval = window.yaxis_tickinterval_saved;
    window.plot.axes.yaxis.numberTicks = window.yaxis_numberticks_saved;
};


function normal_scaleY2(){
    window.plot.axes.y2axis.min = window.y2axis_min_saved;
    window.plot.axes.y2axis.max = window.y2axis_max_saved;
    window.plot.axes.y2axis.tickInterval = window.y2axis_tickinterval_saved;
    window.plot.axes.y2axis.numberTicks = window.y2axis_numberticks_saved;
};


function normalizeY2(){

    if (!window.zoomedin)
        return '';

    var plotData =  window.plot.series[1].data;
    var y=[];
    for (i=0;i<plotData.length;i++) {
        if (plotData[i][0] >= window.plot.axes.xaxis.min && plotData[i][0] <= window.plot.axes.xaxis.max)
            y.push(plotData[i][1]);
    }
    var maxtotal_y = Math.max.apply(null, y);
    var mintotal_y = Math.min.apply(null, y);

    if (maxtotal_y==mintotal_y){
        var w=maxtotal_y;
        mintotal_y=w*0.8;
        maxtotal_y=w*1.2;
    }

    window.plot.axes.y2axis.max=maxtotal_y;
    window.plot.axes.y2axis.min=mintotal_y;


    var delta=(window.plot.axes.y2axis.max-window.plot.axes.y2axis.min)*0.1;
    window.plot.axes.y2axis.min=Math.max(window.plot.axes.y2axis.min-delta,0);
    window.plot.axes.y2axis.max=window.plot.axes.y2axis.max+delta;
    window.plot.axes.y2axis.tickInterval=(window.plot.axes.y2axis.max-window.plot.axes.y2axis.min)/5;
    window.plot.axes.y2axis.numberTicks=6;

    if (window.plot.axes.y2axis.max ==0 && window.plot.axes.y2axis.min ==0){
            window.plot.axes.y2axis.min=0;
            window.plot.axes.y2axis.max=1;
    }


};


function switchAxesToNormal(){

    window.norm = false;

    normal_scaleY();
    if (window.y2_profile != '---') {
        normal_scaleY2();
    }

    window.plot.replot();

};


function switchAxesToRescale(){

    window.norm = true;

    normalizeY();
    if (window.y2_profile != '---'){
        normalizeY2();
    }

    window.plot.replot();

};


function switchOnCutoff(){
    window.showoverlay = true;
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
    window.zoomedin = false;
};


function switchOffCutoff(){
    window.showoverlay = false
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
    window.zoomedin = false;
};


function goToRegionsView(){
    setRegionsViewValues(whichSelected(), '', false)
    window.location = '/regions';
};


function whichSelected(){
    var sel = '';
    $('.selected').each(function() {
        sel = $(this).find('td:first').text();
    });
    return sel
};

function switchToForwardStrand(){
    window.reads = 'forward';
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
    window.zoomedin = false;
};

function switchToReverseStrand(){
    window.reads = 'reverse';
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
    window.zoomedin = false;
};

function switchToBothStrands(){
    window.reads = 'all';
    makePlot();
    makeReferenceBar(window.sequences[window.region]);
    $("#zoombuttons").hide();
    window.zoomedin = false;
};