console.log("hi");
var oTable;

function isDefined(x) {
    var undefined;
    return x !== undefined;
}


function ellipsis(text, n) {
	console.log("in ellipsis function");
    if(text.length>n)
        return text.substring(0,n)+"...";
    else
        return text;
}

var myCallbackFunction = function( nRow, aData, iDisplayIndex ) {
							console.log("in row callback function");
                            // the 1 here indicates which column (counting from 0) to truncate:
                            var $cell=$('td:eq(3)', nRow);
                            $cell.text(ellipsis($cell.text(),10));
                            var $cell=$('td:eq(5)', nRow);
                            $cell.text(ellipsis($cell.text(),10));
                            return nRow;
                    }


$(document).ready(function() {

		console.log("hi from document ready function");
        $('#mytable tbody td:not(:first-child)').editable( function( sValue ) {
                var aPos = oTable.fnGetPosition( this );
                var aData = oTable.fnGetData( aPos[0] );
                aData[ aPos[1] ] = sValue;
                return ellipsis(sValue,20);
        }, { data: function(value, settings) {
                var aPos = oTable.fnGetPosition( this );
                var aData = oTable.fnGetData( aPos[0] );
                return aData[aPos[2]];
              }
        ,"onblur": 'submit' } );

        console.log("is myCallbackFunction defined? " + isDefined(myCallbackFunction));

/*
        oTable = $('#mytable').dataTable( {
                    "fnRowCallback": function( nRow, aData, iDisplayIndex ) {
							console.log("in row callback function");
                            // the 1 here indicates which column (counting from 0) to truncate:
                            var $cell=$('td:eq(4)', nRow);
                            $cell.text(ellipsis($cell.text(),2));
                            return nRow;
                    }
            } );

*/

} );