'use client';

import { useEffect, useState } from 'react';
import dynamic from 'next/dynamic';
import { api } from '@/lib/api';

// Plotlyを動的インポート（SSRを無効化）
const Plot = dynamic(() => import('react-plotly.js'), { ssr: false });

interface PerformanceChartProps {
  articleId?: string;
  metrics?: {
    likes_count?: number[];
    page_views_count?: number[];
    dates?: string[];
  };
}

export default function PerformanceChart({ articleId, metrics: propsMetrics }: PerformanceChartProps) {
  const [chartData, setChartData] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  
  useEffect(() => {
    const loadMetrics = async () => {
      if (!articleId) {
        // プロップスからメトリクスを使用
        updateChartData(propsMetrics);
        return;
      }
      
      setLoading(true);
      try {
        const metricsData = await api.getArticleMetrics(articleId);
        const currentMetrics = metricsData.current_metrics || metricsData.kpi || {};
        
        // 時系列データがない場合は現在のメトリクスを表示
        const metrics = {
          likes_count: [currentMetrics.likes_count || 0],
          page_views_count: [currentMetrics.page_views_count || 0],
          dates: [new Date().toISOString().split('T')[0]]
        };
        
        updateChartData(metrics);
      } catch (err) {
        console.error('メトリクス取得エラー:', err);
        updateChartData(propsMetrics);
      } finally {
        setLoading(false);
      }
    };
    
    loadMetrics();
  }, [articleId, propsMetrics]);
  
  const updateChartData = (metrics?: {
    likes_count?: number[];
    page_views_count?: number[];
    dates?: string[];
  }) => {
    if (metrics && metrics.dates && metrics.dates.length > 0) {
      const data = [
        {
          x: metrics.dates,
          y: metrics.likes_count || [],
          type: 'scatter',
          mode: 'lines+markers',
          name: 'いいね数',
          line: { color: '#3B82F6' }
        },
        {
          x: metrics.dates,
          y: metrics.page_views_count || [],
          type: 'scatter',
          mode: 'lines+markers',
          name: '閲覧数',
          yaxis: 'y2',
          line: { color: '#10B981' }
        }
      ];
      
      const layout = {
        title: 'パフォーマンス推移',
        xaxis: { title: '日付' },
        yaxis: { title: 'いいね数', side: 'left' },
        yaxis2: {
          title: '閲覧数',
          overlaying: 'y',
          side: 'right'
        },
        hovermode: 'x unified',
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { color: '#6B7280' },
        legend: { x: 0, y: 1 }
      };
      
      setChartData({ data, layout });
    } else {
      // デフォルトの空データ
      setChartData({
        data: [{
          x: [],
          y: [],
          type: 'scatter',
          mode: 'lines+markers',
          name: 'いいね数'
        }],
        layout: {
          title: 'パフォーマンス推移',
          xaxis: { title: '日付' },
          yaxis: { title: 'いいね数' }
        }
      });
    }
  };
  
  if (loading || !chartData) {
    return (
      <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6">
        <div className="flex items-center justify-center h-64">
          <p className="text-gray-600 dark:text-gray-400">データを読み込んでいます...</p>
        </div>
      </div>
    );
  }
  
  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6">
      <Plot
        data={chartData.data}
        layout={chartData.layout}
        style={{ width: '100%', height: '400px' }}
        config={{ responsive: true, displayModeBar: false }}
      />
    </div>
  );
}

