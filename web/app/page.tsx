'use client';

import { useEffect, useState } from 'react';
import { api, ArticleListItem } from '@/lib/api';
import ArticleStatusCard from '@/components/ArticleStatusCard';
import ApprovalCard from '@/components/ApprovalCard';
import KPICard from '@/components/KPICard';
import PerformanceChart from '@/components/PerformanceChart';
import LoadingSpinner from '@/components/LoadingSpinner';
import Link from 'next/link';

interface DashboardSummary {
  total_articles: number;
  total_likes: number;
  total_views: number;
  avg_engagement_rate: number;
  pending_approval_count: number;
}

export default function Dashboard() {
  const [articles, setArticles] = useState<ArticleListItem[]>([]);
  const [pendingArticles, setPendingArticles] = useState<ArticleListItem[]>([]);
  const [summary, setSummary] = useState<DashboardSummary | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [filter, setFilter] = useState<'all' | 'published' | 'draft'>('all');
  const [sortBy, setSortBy] = useState<'date' | 'likes' | 'views'>('date');

  useEffect(() => {
    loadDashboardData();
  }, []);

  const loadDashboardData = async () => {
    try {
      setLoading(true);
      setError(null);
      
      // 並列でデータを取得（エラー時は空配列/デフォルト値を設定）
      const [articlesData, pendingData, summaryData] = await Promise.all([
        api.listArticles().catch(() => []),
        api.getPendingApprovalArticles().catch(() => []),
        api.getDashboardSummary().catch(() => ({
          total_articles: 0,
          total_likes: 0,
          total_views: 0,
          avg_engagement_rate: 0.0,
          pending_approval_count: 0
        }))
      ]);
      
      setArticles(articlesData);
      setPendingArticles(pendingData);
      setSummary(summaryData);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'データの読み込みに失敗しました';
      setError(errorMessage);
      console.error('ダッシュボードデータ読み込みエラー:', err);
    } finally {
      setLoading(false);
    }
  };

  // フィルタリングとソート
  const filteredAndSortedArticles = articles
    .filter(article => {
      if (filter === 'published') {
        // 投稿済み（qiitaUrlがある記事は別の方法で判定が必要）
        return true; // 暫定
      } else if (filter === 'draft') {
        return true; // 暫定
      }
      return true;
    })
    .sort((a, b) => {
      if (sortBy === 'date') {
        const dateA = a.updatedAt ? new Date(a.updatedAt).getTime() : 0;
        const dateB = b.updatedAt ? new Date(b.updatedAt).getTime() : 0;
        return dateB - dateA;
      }
      // likes/viewsのソートは将来的に実装
      return 0;
    });

  return (
    <div className="space-y-6">
      {/* ヘッダー */}
      <div className="flex justify-between items-center">
        <h1 className="text-3xl font-bold text-gray-900 dark:text-white">ダッシュボード</h1>
        <Link
          href="/new"
          className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
        >
          新規記事作成
        </Link>
      </div>

      {error && (
        <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          {error}
        </div>
      )}

      {loading ? (
        <LoadingSpinner />
      ) : (
        <>
          {/* KPIサマリー */}
          {summary && (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
              <KPICard
                title="総記事数"
                value={summary.total_articles}
                subtitle="投稿済み記事"
              />
              <KPICard
                title="総いいね数"
                value={summary.total_likes.toLocaleString()}
              />
              <KPICard
                title="総閲覧数"
                value={summary.total_views.toLocaleString()}
              />
              <KPICard
                title="平均エンゲージメント率"
                value={`${summary.avg_engagement_rate}%`}
                subtitle="いいね数/閲覧数"
              />
            </div>
          )}

          {/* 承認待ち記事セクション */}
          {pendingArticles.length > 0 && (
            <div>
              <h2 className="text-xl font-semibold text-gray-900 dark:text-white mb-4">
                承認待ち記事 ({pendingArticles.length})
              </h2>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                {pendingArticles.map((article) => (
                  <ApprovalCard
                    key={article.id}
                    article={{
                      id: article.id,
                      title: article.title,
                      approvalDeadline: (article as any).approvalDeadline,
                      scheduledPublishDate: (article as any).scheduledPublishDate
                    }}
                  />
                ))}
              </div>
            </div>
          )}

          {/* パフォーマンスグラフ */}
          <div>
            <h2 className="text-xl font-semibold text-gray-900 dark:text-white mb-4">
              パフォーマンス推移
            </h2>
            <PerformanceChart />
          </div>

          {/* フィルタリング・ソート */}
          <div className="flex justify-between items-center">
            <div className="flex space-x-4">
              <select
                value={filter}
                onChange={(e) => setFilter(e.target.value as 'all' | 'published' | 'draft')}
                className="px-3 py-2 border border-gray-300 rounded-md dark:bg-gray-700 dark:border-gray-600 dark:text-white"
              >
                <option value="all">すべて</option>
                <option value="published">投稿済み</option>
                <option value="draft">下書き</option>
              </select>
              <select
                value={sortBy}
                onChange={(e) => setSortBy(e.target.value as 'date' | 'likes' | 'views')}
                className="px-3 py-2 border border-gray-300 rounded-md dark:bg-gray-700 dark:border-gray-600 dark:text-white"
              >
                <option value="date">日付順</option>
                <option value="likes">いいね数順</option>
                <option value="views">閲覧数順</option>
              </select>
            </div>
          </div>

          {/* 記事一覧 */}
          <div>
            <h2 className="text-xl font-semibold text-gray-900 dark:text-white mb-4">
              記事一覧 ({filteredAndSortedArticles.length})
            </h2>
            {filteredAndSortedArticles.length === 0 ? (
              <div className="text-center py-12">
                <p className="text-gray-600 dark:text-gray-400 mb-4">記事がありません</p>
                <Link
                  href="/new"
                  className="inline-block px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
                >
                  最初の記事を作成
                </Link>
              </div>
            ) : (
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
                {filteredAndSortedArticles.map((article) => (
                  <ArticleStatusCard key={article.id} article={article} />
                ))}
              </div>
            )}
          </div>
        </>
      )}
    </div>
  );
}
