function construct_prediction_datatable(name, dataTitle, data, sur_cor, high_percentile, low_percentile, stage, cancer, switch_dict, selected_miRNA, DE_info_dict)
{
    //// Dynamically add/remove datatable column&data using ajax.
    // without if condition will raise error about reinitialise(Cannot read property 'aDataSort' of undefined)
    if ( $.fn.DataTable.isDataTable( `#${name}` ) ) {
        $(`#${name}`).DataTable().destroy();
        $(`#${name}`).empty();
    };
    // colspan column
    var headers = '<thead><tr><th colspan="1">Name</th>'
    // <th colspan="2">Title</th><th colspan="2">1</th><th colspan="2">2</th><th colspan="2">3</th><th colspan="2">4</th><th colspan="2">5</th><th colspan="2">6</th><th colspan="2">7</th><th colspan="2">8</th></tr>';
    // default columns
    var columns = [{ title:`${dataTitle}`, data: "name"}];
    var columnDefs = []
    var target_index = 1;
    // create survival analysis columns and columnDefs
    if(switch_dict.survival == true){
        headers += '<th colspan="2">survival analysis</th>';
        if(sur_cor == "None"){
            columns.push({ title:"survival analysis logrank p-value", data: "logrank_p_value"})
        }
        else{
            columns.push({ title:`${sur_cor} correction p-value`, data: `${sur_cor}`})
        }
        columns.push(
            { title:"survival analysis detail"},
        );
        columnDefs.push(
            {
                // 指定第一列，從0開始，0表示第一列，1表示第二列……
                targets: target_index ,
                autoWidth: true,
                render: function(data, type, row, meta) {
                    // console.log(typeof data);
                    if (type === 'display' && typeof data === 'number') {
                        return Number(data.toFixed(6)).toExponential(3);
                    }
                    return data;
                },
            },
            {
                // 指定第一列，從0開始，0表示第一列，1表示第二列……
                targets: target_index + 1,
                autoWidth: true,
                render: function(data, type, row, meta) {
                    // if id has space or : or . it can not be use
                    // var title = meta.settings.aoColumns[meta.col-2].sTitle;
                    // console.log(title)
                    // var name = row[title];
                    var name = row["name"]; // 後面的key是用data:後的名稱
                    return `<a href="/survival_analysis/detail/?cancer=${cancer}&type=${dataTitle}&name=${name}&hp=${high_percentile}&lp=${low_percentile}&stage=${stage}" target="_blank">view</a>`
                },
            },
        );
        target_index += 2;
    }
    // create miRNA columns and columnDefs
    if(switch_dict.miRNA == true){
        headers += `<th colspan="${selected_miRNA.length}">miRNA</th>`;
        selected_miRNA.forEach((miRNA) => {
            columns.push(
                { title:`${miRNA}`, data: `${miRNA}`}
            );
        });
        target_index += selected_miRNA.length;
    }
    // create DE columns and columnDefs
    if(switch_dict.DE == true){
        headers += `<th colspan="4">DE</th>`;
        columns.push(
            { title:`${DE_info_dict["condition2"]} Avg FPKM(cond2)`, data: `avg_s_FPKM`},
            { title:`${DE_info_dict["condition1"]} Avg FPKM(cond1)`, data: `avg_f_FPKM`},
            { title:`Fold Change(cond2/cond1)`, data: `foldchange`},
            { title:`q-value`, data: `q-value`},
        );
        columnDefs.push(
            {   
                // 指定第一列，從0開始，0表示第一列，1表示第二列……
                targets: [target_index, target_index + 1, target_index + 2, target_index + 3],
                autoWidth: true,
                render: function(data, type, row, meta) {
                    // console.log(target_index)
                    // console.log(typeof data);
                    if (type === 'display' && typeof data === 'number') {
                        return Number(data.toFixed(6)).toExponential(3);
                    }
                    return data;
                },
            },
        );
        target_index += 4;
    }
    // $(`#${name}`).append(headers);
    $(`#${name}`).DataTable({
        // 'scrollX':true, will split header and data, it will cause header and data misalign
        // datatable size will auto Resizes with the window
        "initComplete": function (settings, json) {  
            $(`#${name}`).wrap("<div style='overflow:auto; width:100%;position:relative;'></div>");            
        },
        // bSort: false,
        order: [[0, 'asc']],
        destroy : true,
        data: data,
        columns: columns,
        columnDefs: columnDefs,
        dom: 'Bfrtip',
        buttons: [{
            extend: 'csv',
            text: 'Export CSV',
            title: `${cancer}_${dataTitle}_${stage}_${high_percentile}_${low_percentile}`,
            exportOptions: {
                columns: [0, 1]
                // columns: [0]
            }
        }],
        
    });
}

function Tag_len_to_adjust_set() {
    var previousLength = 0;
    // Set up an interval to check the length at regular intervals (every 1000 milliseconds)
    var intervalId = setInterval(function() {
        var currentLength = $('#miRNA_input_area').jsonTagEditor('getTags')[0].tags.length;
        // Check if the length has changed
        if (currentLength !== previousLength) {
            if (currentLength < 2) {
                $("#set_operation_area").css("display", "none"); 
            } else {
                $("#set_operation_area").css("display", "block");
                if(currentLength != 2){
                    $(".diff_set").hide();
                    $("#miRNA_union").prop("checked", true);
                }
                else{
                    $(".diff_set").show();
                }
            }
            // Update the previous length for the next check
            previousLength = currentLength;
        }
    }, 200); // Adjust the interval time as needed (e.g., every second)
}

function showTab(tabId, clickedTab) {
    // Hide all tab contents
    var tabs = document.getElementsByClassName('tab-content');
    for (var i = 0; i < tabs.length; i++) {
      tabs[i].style.display = 'none';
    }

    // Remove 'active' class from all tabs
    var allTabs = document.getElementsByClassName('nav-link');
    for (var i = 0; i < allTabs.length; i++) {
      allTabs[i].classList.remove('active');
    }

    // Show the selected tab content
    document.getElementById(tabId).style.display = 'block';

    // Add 'active' class to the clicked tab
    clickedTab.classList.add('active');
}

function nextTab(tabId, activeTab) {
    var activeId = document.getElementById(activeTab);
    // Hide all tab contents
    var tabs = document.getElementsByClassName('tab-content');
    for (var i = 0; i < tabs.length; i++) {
      tabs[i].style.display = 'none';
    }

    // Remove 'active' class from all tabs
    var allTabs = document.getElementsByClassName('nav-link');
    for (var i = 0; i < allTabs.length; i++) {
      allTabs[i].classList.remove('active');
    }

    // Show the selected tab content
    document.getElementById(tabId).style.display = 'block';

    activeId.classList.add('active');

    activeId.classList.remove('disabled');
}

$(document).ready(function(){
    var dataTitle = 'genes'
    $(".btn-toggle").click(function () {
        // Remove the 'active' class from all buttons
        $(".btn-toggle").removeClass("active");
        $(".btn-toggle").addClass("notActive");

        // Add the 'active' class to the clicked button
        $(this).removeClass("notActive");
        $(this).addClass("active");

        // Update the input value
        dataTitle = $(this).data("title");
        $("#DE_level").val(dataTitle);
    });


    Tag_len_to_adjust_set();

    // for miRNA screener autocomplete data
    var homo_miRNA_list = document.getElementById('homo_miRNA_list').getAttribute('data-json').replace(/[\[\]\'()]/g, '').split(",");
    var miRNA_type = "homo";
    var miRNA_list;
    // miRNA type switch
    switch (miRNA_type) {
        case 'homo':
            miRNA_list = homo_miRNA_list;
            break;
        default:
            alert('沒有符合的條件');
    }

    // miRNA screener setting
    $('#miRNA_input_area').jsonTagEditor({
        autocomplete: {
            minLength: 7,
            delay: 0, // show suggestions immediately
            position: { collision: 'flip' }, // automatic menu position up/down
            source: miRNA_list
        },
        forceLowercase: false,
        placeholder: 'miRNA name(s)'
    });

    $('#search_btn').click(function(){
        $('#output_div').hide()
        $('#user_input_miRNA_screener').hide()
        $('#user_input_survival_analysis').hide()
        $('#user_input_de_screener').hide()
        $('#enrichment_analysis_user_input').hide();
		$('#enrichment_analysis_output').hide();
        
        // get input cancer
        var select_elements = document.querySelectorAll('[name="select_element"]');
        var select_cancer = [];
        for(var i=0; i<select_elements.length; i++) {
            select_cancer.push(select_elements[i].value);
            console.log(i,select_elements[i].value);
        }
        
        //// p-value screener
        var survival_switch = $("#switch_survival").prop('checked'); //true or false
        // get stage
        var stage = document.getElementById("stage_select").value;
        // get percentile
        var high_percent = document.getElementById("high_percentile").value;
        var low_percent = document.getElementById("low_percentile").value;
        // get input pvalue
        var p_value = document.getElementById("p_value").value; 
        // get correction method
        var sur_cor = document.getElementById("c_type_sur").value;
        //// miRNA screener
        var miRNA_switch = $("#switch_miRNA").prop('checked');
        $('#miRNA_input_area').jsonTagEditor('getTags')[0].tags.length;
        var selected_miRNA = [];
        for(var i = 0; i < $('#miRNA_input_area').jsonTagEditor('getTags')[0].tags.length; i++){
            selected_miRNA.push($('#miRNA_input_area').jsonTagEditor('getTags')[0].tags[i].value);
        }
        var miRNA_set = document.querySelector('input[name="miRNA_set"]:checked').value;
        //// DE screener
        var DE_switch = $('#switch_DE').prop('checked');
        var DE_filter_elements = document.querySelectorAll('[name="DE_filter_elements"]');
        var DE_filter = [];
        for(var i=0; i<DE_filter_elements.length; i++) {
            DE_filter.push(DE_filter_elements[i].value);
        }
        var DE_info_dict = new Object();
        DE_info_dict.condition1 = document.getElementById("condition1_pre").value.split('|')[0];
        DE_info_dict.condition2 = document.getElementById("condition2_pre").value.split('|')[0];

        var switch_dict = new Object();
        console.log(survival_switch, miRNA_switch, DE_switch);
        switch_dict.survival = survival_switch;
        switch_dict.miRNA = miRNA_switch;
        switch_dict.DE = DE_switch;
        var switch_string = JSON.stringify(switch_dict)
        
        if( survival_switch==false && miRNA_switch==false && DE_switch==false){
            swal('Please open at least 1 screener');
        }
        else if(select_cancer.includes('') == false ){
            $.ajax({
                headers: { 'X-CSRFToken': csrf_token },
                type: 'POST',
                // url:'/survival_analysis/cal_pvalue/',
                url:'/screener/screener_cal_result_gene/',
                dataType : 'json',
                data : {
                    'switch_dict':switch_string,
                    'type':dataTitle,
                    'cancer':select_cancer[0],
                    // p value
                    'stage': stage,
                    'high_percent': high_percent,
                    'low_percent': low_percent,
                    "pvalue":p_value,
                    "cor_method_sur":sur_cor,
                    // miRNA
                    'selected_miRNA[]':selected_miRNA,
                    'miRNA_set':miRNA_set,
                    // DE
                    'DE_filter[]':DE_filter,
                },
                beforeSend:function(){
                    var count=0
                    tID= setInterval(timedCount , 50);
                        function timedCount() {
                        count=count+0.05;
                        swal({
                            title: "Running...",
                            text: "It may take several minutes.\nPlease be patient.\n \nRunning time: "+parseInt(count)+" seconds\nClick anywhere of the page \nif the running time does not change",                       
                            button: false,
                        });
                    };
                },     
                success:function(response){
                    swal.close();
                    $('#output_div').show()
                    clearInterval(tID);
                    delete tID

                    $(`#output_table`).html();
                    $('#gap').html('<br><br><br>')
                    construct_prediction_datatable(
                        "output_table", dataTitle, response.result, sur_cor,
                        high_percent, low_percent, stage, 
                        select_cancer[0], switch_dict, selected_miRNA ,DE_info_dict);
                    $('#stage_td').text(stage);
                    $('#input_type_td').text(dataTitle);
                    $('#primary_site_td').text(select_cancer[0]);
                    $('#high_percent_td').text(high_percent + "%");
                    $('#low_percent_td').text(low_percent + "%");
                    $('#correction_td').text(sur_cor);
                    $('#input_pvalue_td').text(p_value);
                    $('#output_message').text(`${response.result.length} ${dataTitle} met the input criteria`);
                    $('#screener_type_td').text(response.screener_type);
                    $('#selected_miRNA_td').text(selected_miRNA);
                    $('#miRNA_set_td').text(miRNA_set);
                    $('#condition2_td').text(DE_info_dict.condition1);
                    $('#condition1_td').text(DE_info_dict.condition2);
                    $('#fold_change_td').text(DE_filter[2] + DE_filter[3]);
                    $('#test_td').text(DE_filter[4]);
                    $('#direction_td').text(DE_filter[5]);
                    $('#de_correction_td').text(DE_filter[6]);
                    $('#q-value_td').text(DE_filter[7]);

                    console.log(DE_filter)
                    if(survival_switch == true){
                        $('#user_input_survival_analysis').show()
                    }

                    if(miRNA_switch == true){
                        $('#user_input_miRNA_screener').show()
                    }
                    
                    if(DE_switch == true){
                        $('#user_input_de_screener').show()
                    }

                    nextTab('main_output', 'main_output_nav')
                    $('#enrichment_analysis_nav').removeClass('disabled')
                },
                error:function(xhr, ajaxOptions, thrownError){ 
                    alert(thrownError);
                    clearInterval(tID);
                    delete tID
                    swal.close();
                    swal('error')
                }       
            });    
        }
        else{
            swal('Input Error!')
        }                                                                                                       
    });
    
    //// switch for every screener
    $('#switch_survival').change(function() {
        if ($(this).prop('checked')) {
            $("#survival_analysis_filter").css("display", "block");
        } else {
            $("#survival_analysis_filter").css("display", "none"); 
        }
    });

    $('#switch_miRNA').change(function() {
        if ($(this).prop('checked')) {
            $("#mirna_screener_filter").css("display", "block");
        } else {
            $("#mirna_screener_filter").css("display", "none");
        }
    });

    $('#switch_DE').change(function() {
        if ($(this).prop('checked')) {
            $("#de_screener_filter").css("display", "block");
        } else {
            $("#de_screener_filter").css("display", "none"); 
        }
    });
    
});

