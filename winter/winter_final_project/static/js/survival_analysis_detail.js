function download_csv(){
    window.location.href = `{% static 'cancer_project/static/data/csv_result/Survival_Profile_${primary_site.tostring()}_${High_Percentile.tostring()}_${Low_Percentile.tostring()}.csv' %}`;
}
$(document).ready(function(){
    var project = $("#project").text().trim();
    var primary_site = $("#primary_site").text().trim();
    var input_type = $("#input_type").text().trim();
    var Low_Percentile = $("#Low_Percentile").text().trim();
    var High_Percentile = $("#High_Percentile").text().trim();
    var stage = $("#stage").text().trim();
    var name = $("#name").text().trim();
    var days = $("#days").text().trim();
    var dynamicURL = `../../static/data/csv_result/Survival_Profile_${primary_site}_${name}_${High_Percentile}_${Low_Percentile}.csv`;
    document.getElementById('download_link').href = dynamicURL;

    $('#download_btn').click(function(){

    });

    $('#submit_plot').click(function(){
        project = $("#project").text().trim();
        primary_site = $("#primary_site").text().trim();
        input_type = $("#input_type").text().trim();
        Low_Percentile = $("#low_percent").val();
        High_Percentile = $("#high_percent").val();

        stage = document.getElementById("stage_select").value;
        name = $("#name").text().trim();
        
        var dynamicURL = `../../static/data/csv_result/Survival_Profile_${primary_site}_${name}_${High_Percentile}_${Low_Percentile}.csv`;
        document.getElementById('download_link').href = dynamicURL;

        console.log(stage);
        console.log(days);
        console.log(Low_Percentile);
        console.log(High_Percentile);
        $.ajax({
            headers: { 'X-CSRFToken': csrf_token },
            type: 'POST',
            url:'/survival_analysis/sur_plot/',
            dataType : 'json',
            data : {
                'project':project, 
                'primary_site': primary_site,
                'type':input_type,
                'name':name,
                'stage': stage,
                'High_Percentile': High_Percentile,
                'Low_Percentile': Low_Percentile,
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
                console.log(response.img_str);
                $("#sur_img").attr("src", '');
                $("#sur_img").attr("src", `data:image/png;base64,${response.img_str}`);
                clearInterval(tID);
                delete tID
                console.log(response.img_str);
            },
            error:function(xhr, ajaxOptions, thrownError){ 
                clearInterval(tID);
                delete tID
                swal.close();
                swal('error')
            }       
        });                                                                                                      
    });

    
});