var x = [1,2,3,4,5,6,7,8,9];
var y = [6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25,8.5];
var stdY = [0,0,0,0,0,0,0,0,0];
var xNodeId = ["n_0","n_7","n_15","n_41","n_26","n_54","n_69","n_116","n_88"];
var xLabels = ["E6.5_0","E6.75_0","E7.0_0","E7.25_0","E7.5_0","E7.75_0","E8.0_0","E8.25_0","E8.5_0"];
let w;
let c;

function rotate() { 
  setTimeout(function (){ 
    cy.nodes().positions(
      function( i, node ){ 
        return { x: node.position('y')*3, y: -1*node.position('x')*1.1
      }; 
  }); 
  cy.fit() 
  },100); 
}

(function(){
  document.addEventListener('DOMContentLoaded', function(){
    let $$ = selector => Array.from( document.querySelectorAll( selector ) );
    let $ = selector => document.querySelector( selector );

    let tryPromise = fn => Promise.resolve().then( fn );

    let toJson = obj => obj.json();
    let toText = obj => obj.text();

    let cy;

    let $dataset = $('#data');
    let getDataset = name => fetch(`datasets/${name}`).then( toJson );
    let applyDataset = dataset => {
      // so new eles are offscreen
      cy.zoom(0.001);
      cy.pan({ x: -9999999, y: -9999999 });

      // replace eles
      cy.elements().remove();
      cy.add( dataset );
    }
    let applyDatasetFromSelect = () => Promise.resolve( $dataset.value ).then( getDataset ).then( applyDataset );

    let $color = $('#color');
    let getStylesheet = name => {
      let convert = res => toJson(res);

      return fetch(`stylesheets/base.json`).then( convert );
    };
    let applyStylesheet = stylesheet => {

      cy.style().fromJson( stylesheet ).update();

      if ($color.value == "stage"){
        c = "mapData("+$color.value+",6.5,8.5,blue,red)"
      }else{
        c = "mapData("+$color.value+",0,1,blue,red)"
      }
      cy.style().selector("node").style({'background-color':c}).update()

      plotLayout["title"]["text"] = $color.value
      let ii = 0;
      for (let i = 0; i < 9; i++){
        if (xLabels[i] != "notSet"){
          ii = xNodeId[i];
          y[i] = cy.getElementById(ii).data($color.value)  
        }
      }
      plotData["y"] = y
      Plotly.newPlot('myDiv', plotData, plotLayout);
    };
    let applyStylesheetFromSelect = () => Promise.resolve( $color.value ).then( getStylesheet ).then( applyStylesheet );

    let $weights = $('#weights');
    let getWeights = name => {
      let convert = res => toJson(res);

      return fetch(`stylesheets/base.json`).then( convert );
    };
    let applyWeights = weight => {

      // cy.style().fromJson( weight ).update();

      switch ($weights.value){
        case "weightForward": w = "mapData(weightForward,0,1,0,10)"; break;
        case "weightBackward": w = "mapData(weightBackward,0,1,0,10)"; break;
        case "weightCompensatedForward": w = "mapData(weightCompensatedForward,0,1,0,10)"; break;
      }
      cy.style().selector("edge").style({'width':w}).update()
    };
    let applyWeightsFromSelect = () => Promise.resolve( $weights.value ).then( getWeights ).then( applyWeights );

    let layouts = {
      hierarchical: { // replace with your own layout parameters
        name: 'breadthfirst',
        directed: true,
        padding: 30
      }
    };
    let prevLayout;
    let getLayout = name => Promise.resolve( layouts[ name ] );
    let applyLayout = layout => {
      if( prevLayout ){
        prevLayout.stop();
      }

      let l = prevLayout = cy.makeLayout( layout );
      rotate()

      return l.run().promiseOn('layoutstop');
    }
    let applyLayoutFromSelect = () => Promise.resolve( "hierarchical" ).then( getLayout ).then( applyLayout );

    cy = window.cy = cytoscape({
      container: $('#cy'),

      layout: { // replace with your own layout parameters
          name: 'breadthfirst',
          directed: true,
          padding: 30
      }
    });

    tryPromise( applyDatasetFromSelect ).then( applyStylesheetFromSelect ).then( applyWeightsFromSelect ).then( applyLayoutFromSelect );

    $dataset.addEventListener('change', function(){
      tryPromise( applyDatasetFromSelect ).then( applyLayoutFromSelect );
    });

    $color.addEventListener('change', function(){
      tryPromise( applyStylesheetFromSelect ).then( applyWeightsFromSelect );
    });

    $weights.addEventListener('change', function(){
      tryPromise( applyStylesheetFromSelect ).then( applyWeightsFromSelect );
    });

    let $a = $('#tracing');
    cy.bind('click', 'node', function(node) {
      Promise.resolve($a)
      Promise.resolve($weights)
      // console.log(node.cyTarget.predecessors().edges());
      switch($a.value){
        case "false":

          break;
        case "tracing":
          switch ($weights.value){
            case "weightBackward":
              node.cyTarget.incomers().edges().animate({
                style: {
                  lineColor: "red"
                }
              });
              break; 
            case "weightCompensatedBackward":
              node.cyTarget.incomers().edges().animate({
                style: {
                  lineColor: "red"
                }
              });
              break;          
            case "weightForward":
              node.cyTarget.outgoers().edges().animate({
                style: {
                  lineColor: "red"
                }
              });
              break;
            case "weightCompensatedForward":
              node.cyTarget.outgoers().edges().animate({
                style: {
                  lineColor: "red"
                }
              });
              break;          
          }
        case "expression":
          Promise.resolve($color) 
          y[node.cyTarget.data("stageInt")] = node.cyTarget.data($color.value)
          xLabels[node.cyTarget.data("stageInt")] = node.cyTarget.data("cluster")
          xNodeId[node.cyTarget.data("stageInt")] = node.cyTarget.data("id")
          Plotly.newPlot('myDiv', plotData, plotLayout);
          break;
        }
      let cluster = node.cyTarget.data("cluster");
      let annotation = node.cyTarget.data("annotation");
      let stage = node.cyTarget.data("stage");
      document.getElementById('cluster').innerHTML = "<b>Cluster: </b>"+cluster;
      document.getElementById('annotation').innerHTML = "<b>Annotation: </b>"+annotation;
      document.getElementById('stage').innerHTML = "<b>Stage: </b>E"+stage;
    });
    cy.bind('cxttap', 'node', function(node) {
      // console.log(node.cyTarget.predecessors().edges());
      node.cyTarget.incomers().edges().animate({
        style: {
          lineColor: "lightgrey"
        }
      });
      node.cyTarget.outgoers().edges().animate({
        style: {
          lineColor: "lightgrey"
        }
      });
    });

    // cy.on("click", "node", (evt) => {evt.cyTarget.ancestors().animate({
    //   style: {lineColor: "red"} }) })
    // cy.off("tap", "node", (evt) => {evt.cyTarget.successors().animate({
    //     style: {lineColor: "grey"} }) })
  
    var plotData = [
      {
        x: x,
        y: y,
        error_y: {
          type: 'data',
          array: stdY,
          visible: true
        },
        type: 'scatter'
      }
    ];
    var plotLayout = {
      title:{text: $color.value},
      xaxis:{
        tickmode: "array",
        tickvals: [1,2,3,4,5,6,7,8,9],
        ticktext: xLabels
      }
    }
    Plotly.newPlot('myDiv', plotData, plotLayout);

    let getDataset2 = name => fetch(`datasets/data2.json`).then( toJson );
    let applyDataset2 = dataset => {
      // so new eles are offscreen
      // console.log(dataset)

      var plotUmap = [
        {
          x: dataset["x"],
          y: dataset["y"],
          mode: 'markers',
          type: 'scatter'
        }
      ];
      var plotUmapLayout = {
        title:{text: $color.value},
      }
      Plotly.newPlot('umap', plotUmap, plotUmapLayout);

      return
    }
    function applyDatasetFromSelect2(){ return Promise.resolve( $dataset.value ).then( getDataset2 ).then( applyDataset2 )};
    applyDatasetFromSelect2()


  });
  
})();

// tooltips with jQuery
$(document).ready(() => $('.tooltip').tooltipster());