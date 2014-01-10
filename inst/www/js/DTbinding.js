$(function() {
    $("#select_all_rows").click(function(){
        $(".selectable div table tbody tr").addClass("rowsSelected");
        $(".selectable div table").trigger("change");
    });

    $("#deselect_all_rows").click(function(){
        $(".selectable div table tbody tr").removeClass("rowsSelected");
        $(".selectable div table").trigger("change");
    });
});

$(document).on('click', '.selectable div table tbody tr', function(e){
    //var oTable = $("#DataTables_Table_0").dataTable();
    //var aPos = oTable.fnGetPosition( this );
    //var data = oTable.fnGetData(this);
    //console.log("You clicked on row " + aPos);
    //console.log("data is: " + data);
	var el = $(this);
	if (!e.ctrlKey){
		$(this).siblings().removeClass("rowsSelected");
	}
	$(this).addClass("rowsSelected", this.clicked);
	el.trigger("change");
});	


var selectRowBinding = new Shiny.InputBinding();
$.extend(selectRowBinding, {
	find: function(scope) {
		return $(scope).find(".selectable");
	},
	getValue: function(el){
    tbl = $(el).find("table");
    var out = [];
    $rows = $(tbl).children().children('.rowsSelected');
    if($rows.length == 0) return -1;
    var oTable = $("#DataTables_Table_0").dataTable();

    $rows.each(function(row,v) {
        var aPos = oTable.fnGetPosition( this );
        var data = oTable.fnGetData(this);
        for (var i = 0; i < data.length; i++) {
            if (typeof out[row] === 'undefined') out[row] = [];
            out[row][i] = data[i];
        }
        /*
      $(this).find("td").each(function(cell,v) {
        if (typeof out[row] === 'undefined') out[row] = [];
        out[row][cell] = $(this).text();
      });
*/
    });
    return out;
	},
	setValue: function(el, value) {
	},
	subscribe: function(el, callback) {
		$(el).on("change.selectRowBinding", function(e) {
			callback();
		});
	},
	unsubscribe: function(el) {
	  $(el).off(".selectRowBinding");
	}
});
Shiny.inputBindings.register(selectRowBinding);

