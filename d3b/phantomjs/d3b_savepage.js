
var fs = require('fs');
var system = require('system');
var args = system.args;
var page = require('webpage').create();

page.onError = function (msg, trace) {
    console.log(msg);
    trace.forEach(function(item) {
        console.log('  ', item.file, ':', item.line);
    });
};

function getFileUrl(str) {
  var pathName = fs.absolute(str).replace(/\\/g, '/');
  // Windows drive letter must be prefixed with a slash
  if (pathName[0] !== "/") {
    pathName = "/" + pathName;
  }
  return encodeURI("file://" + pathName);
};

function savedata( filename, data ) {
	fs.write( filename, data, 'w');
}


//var fileUrl = getFileUrl("lastres.log");

page.open( getFileUrl( args[1] ), function(status) {
  console.log("Status: " + status);
  if(status === "success") {
    page.render( args[2] );
	savedata( args[2], page.content );
  }
  phantom.exit();
});


