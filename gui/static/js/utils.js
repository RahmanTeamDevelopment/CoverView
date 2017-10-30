function setRegionsViewValues(region, gene, failed){
    $.ajax({
        url: '/regions',
        contentType: 'application/json;charset=UTF-8',
        data : JSON.stringify({ region: region, gene: gene, failed: failed }),
        type: 'POST',
        success: function(response) {
        },
        error: function(error) {
        }
    });
};


function setProfilesViewValues(region){
    $.ajax({
        url: '/profiles',
        contentType: 'application/json;charset=UTF-8',
        data : JSON.stringify({ region: region }),
        type: 'POST',
        success: function(response) {
        },
        error: function(error) {
        }
    });
};


function setGenesViewValues(gene){
    $.ajax({
        url: '/genes',
        contentType: 'application/json;charset=UTF-8',
        data : JSON.stringify({ gene: gene }),
        type: 'POST',
        success: function(response) {
        },
        error: function(error) {
        }
    });
};

function goToDefaultProfilesView(){
    setProfilesViewValues('');
    window.location = '/profiles';
}

function goToDefaultRegionsView(){
    setRegionsViewValues('', '', false);
    window.location = '/regions';
}

function goToDefaultGenesView(){
    setGenesViewValues('');
    window.location = '/genes';
}